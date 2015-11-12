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

#include "LinearSystem_MultiSplitting.hh"

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

LinearSystem_MultiSplitting::
	LinearSystem_MultiSplitting(Teuchos::RCP<Teuchos::ParameterList> pl, 
	Teuchos::RCP<const CRS_MATRIX>             &A,
	Teuchos::RCP<const MV>                     &b )
{
    
    Teuchos::RCP<Teuchos::ParameterList> multisplit_pl =
        Teuchos::sublist(pl,"MultiSplitting");             

    d_inner_solver = multisplit_pl->get("inner_solver","richardson");
    VALIDATE( d_inner_solver == "richardson" || 
              d_inner_solver =="monte_carlo", 
             "The type of inner solver provided is not valid for Multi-Splitting" );
             
    d_b = b;  
    d_A = A;       

    createPartitions(pl);

}

//---------------------------------------------------------------------------//
/*!
 * \brief Construction of partitions.
 */
//---------------------------------------------------------------------------//
void
LinearSystem_MultiSplitting::createPartitions( Teuchos::RCP<Teuchos::ParameterList> pl )
{    

    Teuchos::RCP<Teuchos::ParameterList> multisplit_pl =
        Teuchos::sublist(pl,"MultiSplitting");

    d_num_blocks   = multisplit_pl->get("num_blocks",10);
    VALIDATE( d_num_blocks >= 2, "Minimal number of partitions is 2" );
 
    d_overlap      = multisplit_pl->get("overlap",0.0);
    VALIDATE( d_overlap>= 0.0 && d_overlap<=1.0, 
            "The percentage of overlapping must be a number between 0 and 1");
                              
    //measure the size of the problem 
    unsigned int N = d_A -> getGlobalNumRows();

    //temporary block size to deremine the entity of the overlapping
    unsigned int size_temp = N / d_num_blocks;
             
    //determine the number of rows that overlaps between adjacent subdomains         
    unsigned int overlapping = d_overlap * size_temp;        
    //std::cout<<"overlapping= "<<overlapping<<std::endl;
    
    //determine the number of rows for each subdomain
    unsigned int block_size = ( N + (d_num_blocks-1)*overlapping ) / d_num_blocks;        
     
    d_partitions.resize(d_num_blocks);

    for(unsigned int i =0; i!=d_num_blocks; ++i)
        d_partitions[i].resize(2);
       
    unsigned int p = 0;
    
    d_partitions[p][0]=0;
    d_partitions[p][1]=block_size - 1;
    p=1;
    
    while(p!=d_num_blocks - 1)
    {
    	d_partitions[p][0] = d_partitions[p-1][1] + 1 - overlapping;
    	d_partitions[p][1] = d_partitions[p-1][1] + 1 - overlapping + block_size;
        p++;
    }
             
    d_partitions[p][0] = d_partitions[p-1][1] + 1 - overlapping;
    d_partitions[p][1] = N-1;
    
    d_partitions_computed = true;

/*    for(unsigned int i =0; i!=d_num_blocks; ++i)
	std::cout<<"partition "<<i<<":\t"<<d_partitions[i][0]<<"\t"<<d_partitions[i][1]<<std::endl;*/

}


//---------------------------------------------------------------------------//
/*!
 * \brief Apply block diagonal scaling to matrix
 */
//---------------------------------------------------------------------------//
Teuchos::RCP<CRS_MATRIX> LinearSystem_MultiSplitting::computeBlockDiagPrec(unsigned int p)
{

    //measure the size of the problem 
    unsigned int N = d_A -> getGlobalNumRows();

    // Create matrix to hold block inverses
    Teuchos::RCP<CRS_MATRIX> invD = Tpetra::createCrsMatrix<SCALAR>(
        d_A->getDomainMap());

    int max_row_elements = d_A->getNodeMaxNumRowEntries();
    unsigned int block_size = d_partitions[p][1] - d_partitions[p][0] + 1;
    
    Teuchos::ArrayRCP<LO>     row_inds(max_row_elements);
    Teuchos::ArrayRCP<SCALAR> row_vals(max_row_elements);
    Teuchos::ArrayRCP<GO>     block_inds(block_size);
    Teuchos::ArrayRCP<SCALAR> block_vals(block_size);
    Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SCALAR> > diag_block(
        new Teuchos::SerialDenseMatrix<LO,SCALAR>(block_size,block_size) );
    Teuchos::SerialDenseSolver<LO,SCALAR> solver;
    Teuchos::RCP<const MAP> indexer = d_A->getColMap();
    size_t num_entries;
    int err;
    for( int index=0; index < d_partitions[p][0]; ++index )
    {
        SCALAR diag_val=0.0;
        d_A->getLocalRowCopy(index,row_inds(),row_vals(),num_entries);
        for( size_t j=0; j<num_entries; ++j )
        {
            if( row_inds[j] == index )
                diag_val = row_vals[j];
        }
        Teuchos::ArrayRCP<GO>     diag_entry_ind(1);
        Teuchos::ArrayRCP<SCALAR> diag_entry_val(1);
        diag_entry_ind[0] = index;
        VALIDATE (diag_val !=0.0, "zero diagonal entry: cannot invert");
        diag_entry_val[0] = 1.0/diag_val;
        invD->insertGlobalValues(index, diag_entry_ind(0,1), diag_entry_val(0,1));
    }
    
    unsigned int block_start = d_partitions[p][0];
    unsigned int block_end = d_partitions[p][1];
    
    // Fill diagonal block
    diag_block->putScalar(0.0);
    for( int i=0; i<block_size; ++i )
    {
        d_A->getLocalRowCopy(block_start+i,row_inds(),row_vals(),num_entries);
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

    for( int index=d_partitions[p][1] + 1; index < N; ++index )
    {
        SCALAR diag_val = 0.0;
        d_A->getLocalRowCopy(index,row_inds(),row_vals(),num_entries);
        for( size_t j=0; j<num_entries; ++j )
        {
            if( row_inds[j] == index )
                diag_val = row_vals[j];
        }
        Teuchos::ArrayRCP<GO>     diag_entry_ind(1);
        Teuchos::ArrayRCP<SCALAR> diag_entry_val(1);
        diag_entry_ind[0] = index;

        VALIDATE (diag_val !=0.0, "zero diagonal entry: cannot invert");
        diag_entry_val[0] = 1.0/diag_val;
        invD->insertGlobalValues(index, diag_entry_ind(0,1), diag_entry_val(0,1));
    }

    invD->fillComplete();
    return invD;
}


splitting
LinearSystem_MultiSplitting::buildSplitting(unsigned int p)
{

    VALIDATE( d_partitions_computed, "partitions are not computed yet");

    splitting split;
    
//    Teuchos::RCP<NODE> node = d-A->getNode();
//    auto A = d_A->clone(node);
    Teuchos::RCP<CRS_MATRIX> DA = Tpetra::createCrsMatrix<SCALAR>(
               d_A->getDomainMap());

    VALIDATE( p<= d_partitions.size(), 
    "Trying to access to a partition with an index bigger than the n. of total partitions." );
	
    unsigned int start = d_partitions[p][0];
    unsigned int end = d_partitions[p][1];
	
    Teuchos::RCP<CRS_MATRIX>  invD = Tpetra::createCrsMatrix<SCALAR>(
    	d_A->getDomainMap());
    invD = computeBlockDiagPrec(p);
	
    Tpetra::MatrixMatrix::Multiply(*invD,false,*d_A,false,*DA,true);

    CHECK( DA->isFillComplete() );

    MV Prb = Tpetra::createCopy(*d_b);
    invD->apply(*d_b,Prb);	

    Teuchos::RCP<MV> Prb_pointer(new MV(Prb));

    split.A = DA;

    Teuchos::RCP<MV> E( new MV(d_A->getDomainMap(),1) );
    Teuchos::ArrayRCP<SCALAR> E_data = E->getDataNonConst(0);
 
    unsigned int N = d_A -> getGlobalNumRows();
    for( unsigned int i=0;  i<N; ++i)
    	E_data[i] = 0.0;

    for( unsigned int i=start;  i<=end; ++i)
        E_data[i] = 1.0;

    if(p!=0)    
    {
		for( unsigned int i=start;  i<=d_partitions[p-1][1]; ++i)
		    E_data[i] = 0.5;
    }   
    
    if(p!=d_num_blocks-1)    
    {
		for( unsigned int i=d_partitions[p+1][0];  i<=end; ++i)
		    E_data[i] = 0.5;
    }   

    split.E = E;
    split.b = Prb_pointer;

    return split;

}


} // namespace alea

