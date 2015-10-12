//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSystemFactory.cc
 * \author Steven Hamilton
 * \brief  Construct EpetraCrsMatrix from ParameterList.
 */
//---------------------------------------------------------------------------//

#include <string>
#include <cmath>
#include <Alea/config.h>

#include "LinearSystemFactory.hh"

// Trilinos includes
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "TpetraExt_MatrixMatrix.hpp"

// Profugus
#include "spn/Eigenvalue_Solver.hh"
#include "spn/Fixed_Source_Solver.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "driver/Manager.hh"
#include "harness/DBC.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Build matrix from ParameterList.
 *
 *  \param pl ParameterList containing information about matrix construction.
 *
 *  Matrix construction is controlled by the following PL entries on the
 *  nested "Problem" sublist:
 *   - matrix_type(string): laplacian, convection-diffusion, matrix_market,
 *                          profugus (no default)
 *   - matrix_size(int):    N>0 (25) (only used if "matrix_type"="laplacian"
 *                          or "convection-diffusion")
 *   - matrix_file(string): valid matrix market filename
 *   - viscosity(SCALAR):   value for viscosity
 *                          (if "matrix_type"="convection-diffusion")
 *   - laplacian_shift(SCALAR): diagonal shift to apply to Laplacian matrix
 *                              (only applies if "matrix_type"="laplacian")
 */
//---------------------------------------------------------------------------//
Teuchos::RCP<LinearSystem>
LinearSystemFactory::buildLinearSystem( Teuchos::RCP<Teuchos::ParameterList> pl )
{
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(pl,"Problem");

    VALIDATE( mat_pl->isType<std::string>("matrix_type"),
            "Must specify matrix_type." );
    std::string matrix_type = mat_pl->get<std::string>("matrix_type");
    VALIDATE(matrix_type=="laplacian"            ||
             matrix_type=="convection-diffusion" ||
             matrix_type=="matrix_market"        ||
             matrix_type=="profugus",
             "Invalid matrix_type specified.");

    Teuchos::RCP<CRS_MATRIX> A = Teuchos::null;
    Teuchos::RCP<MV>     b = Teuchos::null;
    if( matrix_type == "laplacian" )
    {
        buildDiffusionSystem(mat_pl,A,b);
    }
    else if( matrix_type == "convection-diffusion" )
    {
        buildConvectionDiffusionSystem(mat_pl,A,b);
    }
    else if( matrix_type == "matrix_market" )
    {
        buildMatrixMarketSystem(mat_pl,A,b);
    }
    else if( matrix_type == "profugus" )
    {
        buildProfugusSystem(mat_pl,A,b);
    }

    A = applyShift(A,mat_pl);
    A = applyScaling(A,b,mat_pl);

    return Teuchos::rcp( new LinearSystem(A,b) );
}

//---------------------------------------------------------------------------//
// Build matrix and RHS for diffusion problem
//---------------------------------------------------------------------------//
void LinearSystemFactory::buildDiffusionSystem(
    Teuchos::RCP<Teuchos::ParameterList> pl,
    Teuchos::RCP<CRS_MATRIX>                 &A,
    Teuchos::RCP<MV>                     &b )
{
    GO N = pl->get("matrix_size",25);

    A = buildLaplacianMatrix(N,pl);

    // Sinusoidal source with homogeneous Dirichlet boundaries
    b = Teuchos::rcp( new MV(A->getDomainMap(),1) );
    const double pi = 4.0*std::atan(1.0);
    Teuchos::ArrayRCP<double> b_data = b->getDataNonConst(0);
    LO Nlocal = b_data.size();
    for( LO i=0; i<Nlocal; ++i )
    {
        double gid = static_cast<double>(b->getMap()->getGlobalElement(i));
        double val = std::sin(pi*gid / static_cast<double>(N-1));
        b->replaceLocalValue(i,0,val);
    }
}

//---------------------------------------------------------------------------//
// Build matrix and RHS for diffusion problem
//---------------------------------------------------------------------------//
void LinearSystemFactory::buildConvectionDiffusionSystem(
    Teuchos::RCP<Teuchos::ParameterList> pl,
    Teuchos::RCP<CRS_MATRIX>                 &A,
    Teuchos::RCP<MV>                     &b )
{
    GO N = pl->get("matrix_size",25);
    SCALAR nu = pl->get("viscosity",1.0);
    Teuchos::RCP<CRS_MATRIX> Ad = buildLaplacianMatrix(N,pl);
    Teuchos::RCP<CRS_MATRIX> Ac = buildConvectionMatrix(N,pl);
    A = Teuchos::rcp_dynamic_cast<CRS_MATRIX>(
        Ad->add(1.0,*Ac,nu,Ad->getDomainMap(),Ad->getRangeMap(),pl),true);
    CHECK( A != Teuchos::null );

    // Sinusoidal source with homogeneous Dirichlet boundaries
    b = Teuchos::rcp( new MV(A->getDomainMap(),1) );
    const double pi = 4.0*std::atan(1.0);
    Teuchos::ArrayRCP<double> b_data = b->getDataNonConst(0);
    LO Nlocal = b_data.size();
    for( LO i=0; i<Nlocal; ++i )
    {
        double gid = static_cast<double>(b->getMap()->getGlobalElement(i));
        double val = std::sin(2*pi*gid/static_cast<double>(N-1));
        b->replaceLocalValue(i,0,std::sin(val));
    }
}

//---------------------------------------------------------------------------//
// Build cell-based 1D finite difference Laplacian matrix on [0,1]
//---------------------------------------------------------------------------//
Teuchos::RCP<CRS_MATRIX>
LinearSystemFactory::buildLaplacianMatrix(
    int N, Teuchos::RCP<Teuchos::ParameterList> pl )
{
    SCALAR h = 1.0 / static_cast<SCALAR>(N);

    SCALAR laplacian_shift = pl->get("laplacian_shift",0.0);
    VALIDATE( laplacian_shift >= 0.0, "laplacian shift must be non-negative");

    // Build comm and map
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    LO base=0;
    Teuchos::RCP<MAP> map( new MAP(N,base,comm) );

    LO numPerRow = 3;
    Teuchos::RCP<CRS_MATRIX> A =
        Teuchos::rcp( new CRS_MATRIX( map, numPerRow ) );

    LO numMyRows = map->getNodeNumElements();

    Teuchos::ArrayRCP<SCALAR> val(3);
    Teuchos::ArrayRCP<GO>     ind(3);
    for( LO i=0; i<numMyRows; ++i )
    {
        LO numThisRow;

        GO gid = map->getGlobalElement( i );
        if( gid == 0 )
        {
            // Dirichlet on left
            ind[0] = 0;
            val[0] = (3.0+laplacian_shift)/(h*h);
            ind[1] = 1;
            val[1] = -1.0/(h*h);
            numThisRow = 2;
        }
        else if( gid == N-1 )
        {
            // Dirichlet on right
            ind[0] = N-2;
            val[0] = -1.0/(h*h);
            ind[1] = N-1;
            val[1] = (3.0+laplacian_shift)/(h*h);
            numThisRow = 2;
        }
        else
        {
            ind[0] = gid-1;
            val[0] = -1.0/(h*h);
            ind[1] = gid;
            val[1] = (2.0+laplacian_shift)/(h*h);
            ind[2] = gid+1;
            val[2] = -1.0/(h*h);
            numThisRow = 3;
        }

        Teuchos::ArrayView<SCALAR> val_view = val(0,numThisRow);
        Teuchos::ArrayView<GO>     ind_view = ind(0,numThisRow);
        A->insertGlobalValues( gid, ind_view, val_view );
    }
    A->fillComplete();

    return A;
}

//---------------------------------------------------------------------------//
// Build cell-based 1D finite difference convection matrix on [0,1]
//---------------------------------------------------------------------------//
Teuchos::RCP<CRS_MATRIX>
LinearSystemFactory::buildConvectionMatrix(
    int N, Teuchos::RCP<Teuchos::ParameterList> pl )
{
    SCALAR h = 1.0 / static_cast<SCALAR>(N);

    // Build comm and map
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    LO base=0;
    Teuchos::RCP<MAP> map( new MAP(N,base,comm) );

    LO numPerRow = 3;
    Teuchos::RCP<CRS_MATRIX> A =
        Teuchos::rcp( new CRS_MATRIX( map, numPerRow ) );

    LO numMyRows = map->getNodeNumElements();

    Teuchos::ArrayRCP<SCALAR> val(2);
    Teuchos::ArrayRCP<GO>     ind(2);
    for( LO i=0; i<numMyRows; ++i )
    {
        LO numThisRow;

        GO gid = map->getGlobalElement( i );
        if( gid == 0 )
        {
            // Dirichlet on left
            ind[0] = 0;
            val[0] = 0.5/h;
            ind[1] = 1;
            val[1] = 0.5/h;
            numThisRow = 2;
        }
        else if( gid == N-1 )
        {
            // Dirichlet on right
            ind[0] = N-2;
            val[0] = 0.5/h;
            ind[1] = N-1;
            val[1] = 0.5/h;
            numThisRow = 2;
        }
        else
        {
            ind[0] = gid-1;
            val[0] = -0.5/h;
            ind[1] = gid+1;
            val[1] = 0.5/h;
            numThisRow = 2;
        }

        Teuchos::ArrayView<SCALAR> val_view = val(0,numThisRow);
        Teuchos::ArrayView<GO>     ind_view = ind(0,numThisRow);
        A->insertGlobalValues( gid, ind_view, val_view );
    }
    A->fillComplete();

    return A;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read matrix from Matrix Market file.
 */
//---------------------------------------------------------------------------//
void LinearSystemFactory::buildMatrixMarketSystem(
        Teuchos::RCP<Teuchos::ParameterList> pl,
        Teuchos::RCP<CRS_MATRIX>                 &A,
        Teuchos::RCP<MV>                     &b )
{
    VALIDATE( pl->isType<std::string>("matrix_filename"),
            "Must specify matrix_filename to build matrix market system.");

    std::string filename = pl->get<std::string>("matrix_filename");

    Teuchos::RCP<NODE> node = KokkosClassic::Details::getNode<NODE>();
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();            
    A = Tpetra::MatrixMarket::Reader<CRS_MATRIX>::readSparseFile(
        filename,comm,node);
	
    Teuchos::RCP<const MAP> map = A->getDomainMap();
        
    if( pl->isType<std::string>("rhs_filename") )
    {
        std::string rhs_file = pl->get<std::string>("rhs_filename");
        b = Tpetra::MatrixMarket::Reader<CRS_MATRIX>::readVectorFile(
            rhs_file,comm,node,map);
    }
        
    else
    {
        // Just make RHS constant
        b = Teuchos::rcp( new MV(map,1) );
        b->putScalar(1.0);
    }

    // If an initial guess is specified, compute the initial residual and
    // use that as the RHS
    if( pl->isType<std::string>("init_guess_filename") )
    {
        std::string x0_file = pl->get<std::string>("init_guess_filename");
        Teuchos::RCP<MV> x0 =
            Tpetra::MatrixMarket::Reader<CRS_MATRIX>::readVectorFile(
                x0_file,comm,node,map);
        A->apply(*x0,*b,Teuchos::NO_TRANS,-1.0,1.0);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read matrix from Matrix Market file.
 */
//---------------------------------------------------------------------------//
void LinearSystemFactory::buildProfugusSystem(
        Teuchos::RCP<Teuchos::ParameterList> pl,
        Teuchos::RCP<CRS_MATRIX>                 &A,
        Teuchos::RCP<MV>                     &b )
{
    // Get name of profugus SPN input file
    VALIDATE( pl->isType<std::string>("profugus_input"),
            "Must specify profugus_input to build Profugus system.");
    std::string profugus_file = pl->get<std::string>("profugus_input");

    // Build manager
    spn::Manager spn_manager;
    spn_manager.setup(profugus_file);

    Teuchos::RCP<const profugus::Solver_Base> solver =
        spn_manager.get_solver();

    // Try to figure out what type of solver we have
    using profugus::TpetraTypes;
    using profugus::Eigenvalue_Solver;
    Teuchos::RCP<const Eigenvalue_Solver<TpetraTypes> > eig_solver =
        Teuchos::rcp_dynamic_cast<const Eigenvalue_Solver<TpetraTypes> >(solver);
    if( eig_solver != Teuchos::null )
    {
        // Get LHS matrix
        const profugus::Linear_System<TpetraTypes> &system =
            eig_solver->get_linear_system();
        A = system.get_Matrix();
        CHECK( A != Teuchos::null );

        // Eigenvalue problem doesn't have an rhs, so build one based on
        //  B*x, where x = [1,1,1,...]^T
        b = system.get_RHS();
        CHECK( b != Teuchos::null );
        b->putScalar(1.0);
        Teuchos::RCP<OP> B = system.get_fission_matrix();
        CHECK( B != Teuchos::null );
        B->apply(*b,*b);
    }
    else
    {
        using profugus::Fixed_Source_Solver;
        Teuchos::RCP<const Fixed_Source_Solver<TpetraTypes> > fixed_solver =
            Teuchos::rcp_dynamic_cast<const Fixed_Source_Solver<TpetraTypes> >(solver);
        if( fixed_solver != Teuchos::null )
        {
            // Get LHS matrix
            const profugus::Linear_System<TpetraTypes> &system =
                fixed_solver->get_linear_system();
            A = system.get_Matrix();
            CHECK( A != Teuchos::null );

            // Get RHS vector
            b = system.get_RHS();
            CHECK( b != Teuchos::null );
        }
        else
        {
            INSIST(false,"Unable to determine profugus solver type.  "
                "Is trilinos_implementation=tpetra?");
        }
    }
}

//---------------------------------------------------------------------------//
/*!
* \brief Apply specified shift to matrix
  */
//--------------------------------------------------------------------------
Teuchos::RCP<CRS_MATRIX> LinearSystemFactory::applyShift(
    Teuchos::RCP<CRS_MATRIX>                 A,
    Teuchos::RCP<Teuchos::ParameterList> pl )
{
    VECTOR diag(A->getDomainMap());
    A->getLocalDiagCopy(diag);

    SCALAR shift = pl->get<SCALAR>("diagonal_shift", 0.0);
     
    Teuchos::ArrayRCP<SCALAR> diag_vec = diag.getDataNonConst();
      
    int max_entries = A -> getNodeMaxNumRowEntries();
    Teuchos::ArrayRCP<LO> inds( max_entries );
    Teuchos::ArrayRCP<SCALAR> vals( max_entries );
    
    A -> resumeFill();
    for ( int i=0; i<diag_vec.size(); ++i )
    {
    	size_t num_entries;
    	A -> getLocalRowCopy( i, inds(), vals(), num_entries );
    	int index = std::find( inds.begin(), inds.end(), i ) - inds.begin();
    	vals[index] = vals[index] + shift * vals[index];
    	Teuchos::ArrayView<SCALAR> vals_view=vals( index, 1 );
    	Teuchos::ArrayView<LO> cols_view=inds( index, 1 );
    	//std::cout<<"size of inds = "<<inds.size()<<" and size of cols = " <<vals.size()<<std::endl;
    	auto changes = A -> replaceLocalValues( i, cols_view, vals_view );
    	CHECK( changes==1 );
    	//std::cout<<num_entries<<std::endl;
     }
        	    	    	            
     A -> fillComplete();
    	    	    	    	                    
     CHECK( A -> isFillComplete() );    
     return A;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply specified scaling to matrix
 */
//---------------------------------------------------------------------------//
Teuchos::RCP<CRS_MATRIX> LinearSystemFactory::applyScaling(
    Teuchos::RCP<CRS_MATRIX>                 A,
    Teuchos::RCP<MV>                     b,
    Teuchos::RCP<Teuchos::ParameterList> pl )
{
    VECTOR diag(A->getDomainMap());
    A->getLocalDiagCopy(diag);

    std::string scale_type = pl->get<std::string>("scaling_type","diagonal");
    VALIDATE( scale_type=="diagonal" ||
              scale_type=="block"    ||
              scale_type=="row"      ||
              scale_type=="column"   ||
              scale_type=="sign"     ||
              scale_type=="file"     ||
              scale_type=="none",
              "scale_type must be one of: none, diagonal, block, row, sign.");

    if( scale_type == "diagonal" )
    {
        Teuchos::ArrayRCP<SCALAR> diag_data = diag.getDataNonConst();
        for( LO i=0; i<diag_data.size(); ++i )
        {
            if( diag_data[i] != 0.0 )
            {
                diag_data[i] = 1.0 / diag_data[i];
            }
            else
            {
                diag_data[i] = 1.0;
            }
        }
        A->leftScale(diag);
        Teuchos::ArrayRCP<SCALAR> b_data = b->getDataNonConst(0);
        for( LO i=0; i<b_data.size(); ++i )
        {
            b_data[i] = b_data[i] * diag_data[i];
            CHECK( !SCALAR_TRAITS::isnaninf(b_data[i]) );
        }
    }
    else if( scale_type == "block" )
    {
        A = applyBlockDiagScaling(A,b,pl);
    }
    else if( scale_type == "row" )
    {
        VECTOR inv_row_sums(A->getDomainMap());
        Teuchos::ArrayRCP<SCALAR> row_sum_data = inv_row_sums.getDataNonConst();

        int max_row_size = A->getNodeMaxNumRowEntries();
        Teuchos::ArrayRCP<SCALAR> row_copy(max_row_size);
        Teuchos::ArrayRCP<LO>     ind_copy(max_row_size);
        size_t num_entries;
        for( LO irow=0; irow<A->getNodeNumRows(); ++irow )
        {
            // Add absolute value of row entries
            A->getLocalRowCopy(irow,ind_copy(),row_copy(),num_entries);
            for( size_t i=0; i<num_entries; ++i )
            {
                row_sum_data[irow] += SCALAR_TRAITS::magnitude(row_copy[i]);
            }

            // Invert row sums
            if( row_sum_data[irow] > 0.0 )
            {
                row_sum_data[irow] = 1.0 / row_sum_data[irow];
            }
            else
            {
                row_sum_data[irow] = 1.0;
            }
        }

        // Modify scaling by sign of diagonal
        Teuchos::ArrayRCP<const SCALAR> diag_data = diag.getData();
        for( int i=0; i<row_sum_data.size(); ++i )
        {
            if( diag_data[i] < 0.0 )
            {
                row_sum_data[i] *= -1.0;
            }
        }

        A->leftScale(inv_row_sums);
        for( size_t i=0; i<b->getLocalLength(); ++i )
        {
            b->replaceLocalValue(i,0,b->getData(0)[i]*row_sum_data[i]);
        }
    }
    else if( scale_type == "sign" )
    {
        // Scale each row by sign of diagonal
        VECTOR diag_sign(A->getDomainMap());
        for( size_t i=0; i<diag_sign.getLocalLength(); ++i )
        {
            diag_sign.replaceLocalValue(i,diag.getData(0)[i] > 0.0 ? 1.0 : -1.0);
        }

        A->leftScale(diag_sign);
        for( size_t i=0; i<b->getLocalLength(); ++i )
        {
            b->replaceLocalValue(i,0,b->getData(0)[i]*diag_sign.getData()[i]);
        }
    }
    
    else if (scale_type == "file")
    {
    	std::string pos_scale = pl->get<std::string>("position_scaling", "left");
    	VALIDATE( (pos_scale == "left" || pos_scale=="right"), "invalid position specified to apply the preconditioner");
    	
    	std::string precond_file = pl->get<std::string>("preconditioner_file","none");  	
    	VALIDATE(precond_file != "none", 
    	"When the scaling_type is 'file' the preconditioner file must be specified ");

	Teuchos::RCP<NODE> node = KokkosClassic::Details::getNode<NODE>();
        Teuchos::RCP<const Teuchos::Comm<int> > comm =
                Teuchos::DefaultComm<int>::getComm();            
    	
    	std::cout<<"Loading preconditioner file"<<std::endl;
	Teuchos::RCP<CRS_MATRIX> Pr = Tpetra::MatrixMarket::Reader<CRS_MATRIX>::readSparseFile(
        precond_file,comm,node);
        
        VALIDATE( (  A->getNodeNumRows() == Pr->getNodeNumRows() ), 
        "System matrix and preconditioner do not have the same number of rows" );
        
        VALIDATE( ( A->getNodeNumCols() == Pr->getNodeNumCols() ), 
        "System matrix and preconditioner do not have the same number of columns" );    
    
    	if ( pos_scale == "left" )
        {
		VALIDATE( ( A->getNodeNumCols() == Pr->getNodeNumRows() ), 
		"System matrix and preconditioner must be compatible for the matrix-matrix multiply" );        
	    
	    	Teuchos::RCP<CRS_MATRIX> C = Tpetra::createCrsMatrix<SCALAR> ( A->getDomainMap());
	    	
	    	//Teuchos::RCP<CRS_MATRIX> C;
	    	//size_t max_nnz = d_A->getNodeMaxNumRowEntries();
	    	//C = Teuchos::rcp( new CRS_MATRIX(d_A->getDomainMap(),max_nnz) );
	    	
	    	Tpetra::MatrixMatrix::Multiply(*Pr, Teuchos::NO_TRANS ,*A, Teuchos::NO_TRANS ,*C, true);
	    	CHECK( C->isFillComplete() );
	    
	    	//Tpetra::MatrixMatrix::Add(C, Teuchos::NO_TRANS , 1.0, Pr, Teuchos::NO_TRANS , 0.0, A);
	    	A=C;
	    	CHECK( A->isFillComplete() );
	    	
	    	MV Prb = Tpetra::createCopy(*b);
	    	
	    	//Pr.multiply(*b, Prb, 1.0, 0.0);
	    	Pr->apply(Prb,*b);
	    	
	    	std::cout<<"Preconditioner applied to matrix and right hand side for left scaling"<<std::endl;
    	}
    	else if ( pos_scale == "right" ) 
    	{
    		VALIDATE( ( A->getNodeNumRows() == Pr->getNodeNumCols() ), 
    		"System matrix and preconditioner must be compatible for the matrix-matrix multiply" ); 
	    	Teuchos::RCP<CRS_MATRIX> C = Tpetra::createCrsMatrix<SCALAR> ( A->getDomainMap());
	    	
	    	//Teuchos::RCP<CRS_MATRIX> C;
	    	//size_t max_nnz = d_A->getNodeMaxNumRowEntries();
	    	//C = Teuchos::rcp( new CRS_MATRIX(d_A->getDomainMap(),max_nnz) );
	    	
	    	Tpetra::MatrixMatrix::Multiply(*A, Teuchos::NO_TRANS ,*Pr, Teuchos::NO_TRANS ,*C, true);
	    	CHECK( C->isFillComplete() );
	    
	    	//Tpetra::MatrixMatrix::Add(C, Teuchos::NO_TRANS , 1.0, Pr, Teuchos::NO_TRANS , 0.0, A);
	    	A=C;
	    	CHECK( A->isFillComplete() );
    		std::cout<<"Preconditioner applied to matrix for right scaling"<<std::endl;
    	}
    }
    
    return A;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply block diagonal scaling to matrix
 */
//---------------------------------------------------------------------------//
Teuchos::RCP<CRS_MATRIX> LinearSystemFactory::applyBlockDiagScaling(
    Teuchos::RCP<CRS_MATRIX>                 A,
    Teuchos::RCP<MV>                     b,
    Teuchos::RCP<Teuchos::ParameterList> pl )
{
    int block_size = pl->get("block_size",1);

    REQUIRE( A->getNodeNumRows() % block_size == 0 );
    LO num_blocks = A->getNodeNumRows() / block_size;

    // Create matrix to hold block inverses
    Teuchos::RCP<CRS_MATRIX> invD = Tpetra::createCrsMatrix<SCALAR>(
        A->getDomainMap());

    int max_row_elements = A->getNodeMaxNumRowEntries();
    Teuchos::ArrayRCP<LO>     row_inds(max_row_elements);
    Teuchos::ArrayRCP<SCALAR> row_vals(max_row_elements);
    Teuchos::ArrayRCP<GO>     block_inds(block_size);
    Teuchos::ArrayRCP<SCALAR> block_vals(block_size);
    Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SCALAR> > diag_block(
        new Teuchos::SerialDenseMatrix<LO,SCALAR>(block_size,block_size) );
    Teuchos::SerialDenseSolver<LO,SCALAR> solver;
    Teuchos::RCP<const MAP> indexer = A->getColMap();
    size_t num_entries;
    int err;
    for( int iblock=0; iblock<num_blocks; ++iblock )
    {
        // Fill diagonal block
        diag_block->putScalar(0.0);
        for( int i=0; i<block_size; ++i )
        {
            int block_start = iblock*block_size;
            int block_end = block_start+block_size;
            A->getLocalRowCopy(block_start+i,row_inds(),row_vals(),num_entries);
            for( size_t j=0; j<num_entries; ++j )
            {
                if( row_inds[j] >= block_start && row_inds[j] < block_end )
                {
                    (*diag_block)(i,row_inds[j]-block_start) = row_vals[j];
                }
            }
        }

        // Invert block
        solver.setMatrix(diag_block);
        err = solver.invert();
        CHECK( 0 == err );

        // Set values in invD
        for( int i=0; i<block_size; ++i )
        {
            int block_start = iblock*block_size;
            int global_i = indexer->getGlobalElement(block_start+i);
            int num_to_insert=0;
            for( int j=0; j<block_size; ++j )
            {
                if( (*diag_block)(i,j) != 0.0 )
                {
                    int global_j = indexer->getGlobalElement(block_start+j);
                    block_inds[num_to_insert] = global_j;
                    block_vals[num_to_insert] = (*diag_block)(i,j);
                    num_to_insert++;
                }
            }
            invD->insertGlobalValues(global_i,block_inds(0,num_to_insert),
                block_vals(0,num_to_insert));
        }
    }
    invD->fillComplete();

    // Form preconditioned matrix
    Teuchos::RCP<CRS_MATRIX> DA = Tpetra::createCrsMatrix<SCALAR>(
        A->getDomainMap());
    Tpetra::MatrixMatrix::Multiply(*invD,false,*A,false,*DA,true);
    CHECK( DA->isFillComplete() );

    MV Prb = Tpetra::createCopy(*b);

    invD->apply(Prb,*b);
	    	 
    return DA;
}

} // namespace alea

