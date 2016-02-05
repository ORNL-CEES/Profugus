//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSystemFactory.cc
 * \author Massimiliano Lupo Pasini
 * \brief  Construct EpetraCrsMatrix from ParameterList.
 */
//---------------------------------------------------------------------------//

#include <string>
#include <cmath>
#include <Alea/config.h>

#include "EigenSystemFactory.hh"

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
#include "spn_driver/Manager.hh"
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
Teuchos::RCP<EigenSystem>
EigenSystemFactory::buildEigenSystem( Teuchos::RCP<Teuchos::ParameterList> pl )
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

    return Teuchos::rcp( new EigenSystem(A,b) );
}

//---------------------------------------------------------------------------//
// Build matrix and RHS for diffusion problem
//---------------------------------------------------------------------------//
void EigenSystemFactory::buildDiffusionSystem(
    Teuchos::RCP<Teuchos::ParameterList> pl,
    Teuchos::RCP<CRS_MATRIX>                 &A,
    Teuchos::RCP<MV>                     &b )
{
    GO N = pl->get("matrix_size",25);

    A = buildLaplacianMatrix(N,pl);

    // Random initial guess
    b = Teuchos::rcp( new MV(A->getDomainMap(),1) );
    Teuchos::ArrayRCP<double> b_data = b->getDataNonConst(0);
    LO Nlocal = b_data.size();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    for( LO i=0; i<Nlocal; ++i )
    {
        double val = dis(gen);
        b->replaceLocalValue(i,0,val);
    }
}

//---------------------------------------------------------------------------//
// Build matrix and RHS for diffusion problem
//---------------------------------------------------------------------------//
void EigenSystemFactory::buildConvectionDiffusionSystem(
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

    // Random initial guess
    b = Teuchos::rcp( new MV(A->getDomainMap(),1) );
    Teuchos::ArrayRCP<double> b_data = b->getDataNonConst(0);
    LO Nlocal = b_data.size();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    for( LO i=0; i<Nlocal; ++i )
    {
        double val = dis(gen);
        b->replaceLocalValue(i,0,val);
    }
}

//---------------------------------------------------------------------------//
// Build cell-based 1D finite difference Laplacian matrix on [0,1]
//---------------------------------------------------------------------------//
Teuchos::RCP<CRS_MATRIX>
EigenSystemFactory::buildLaplacianMatrix(
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
EigenSystemFactory::buildConvectionMatrix(
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
void EigenSystemFactory::buildMatrixMarketSystem(
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
        
    if( pl->isType<std::string>("initial_guess") )
    {
        std::string rhs_file = pl->get<std::string>("initial_guess");
        b = Tpetra::MatrixMarket::Reader<CRS_MATRIX>::readVectorFile(
            rhs_file,comm,node,map);
    }
        
    else
    {
    	// Random initial guess
    	b = Teuchos::rcp( new MV(A->getDomainMap(),1) );
    	Teuchos::ArrayRCP<double> b_data = b->getDataNonConst(0);
    	LO Nlocal = b_data.size();
    	std::random_device rd;
    	std::mt19937 gen(rd());
    	std::uniform_real_distribution<> dis(0, 1);

    	for( LO i=0; i<Nlocal; ++i )
    	{
        	double val = dis(gen);
        	b->replaceLocalValue(i,0,val);
    	}
    }

}

//---------------------------------------------------------------------------//
/*!
 * \brief Read matrix from Matrix Market file.
 */
//---------------------------------------------------------------------------//
void EigenSystemFactory::buildProfugusSystem(
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
Teuchos::RCP<CRS_MATRIX> EigenSystemFactory::applyShift(
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

} // namespace alea

