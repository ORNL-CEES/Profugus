//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSystemFactory_MultiSplitting.hh
 * \author Massimiliano Lupo Pasini
 * \brief  Construct MultiSplitting.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_LinearSystem_MultiSplitting_hh
#define mc_solvers_LinearSystem_MultiSplitting_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "LinearSystem.hh"
#include "AleaTypedefs.hh"
#include "MultiSplitting_Utils.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class LinearSystemFactory
 * \brief Load/construct matrix.
 *
 * This class provides a simple interface to build a linear system.  Several
 * simple matrices can be constructed from the provided ParameterList or
 * a matrix can be loaded from a Matrix Market file.  A right hand side for
 * the linear system will also be constructed, either by generating a physical
 * source term where appropriate or by assigning a dummy vector if no
 * meaningful vector is available.
 */
//---------------------------------------------------------------------------//
class LinearSystem_MultiSplitting
{
 
  public:

    LinearSystem_MultiSplitting(Teuchos::RCP<Teuchos::ParameterList>);

    void createPartitions(Teuchos::RCP<Teuchos::ParameterList>);    

    splitting buildSplitting(Teuchos::RCP<Teuchos::ParameterList>, unsigned int);
       
    void buildSystem(Teuchos::RCP<Teuchos::ParameterList>,   
        Teuchos::RCP<CRS_MATRIX> &,
        Teuchos::RCP<MV>         &);
        
    inline std::string getInnerSolverType(){ return d_inner_solver; };

  private:

    //Data
    Teuchos::RCP<CRS_MATRIX> d_A;
    Teuchos::RCP<MV> d_b; 
    
    unsigned int d_num_blocks;
    SCALAR d_overlap;
    std::string d_inner_solver;  
    
    Teuchos::ArrayRCP< ENDPOINTS > d_partitions;


    //Methods
	void buildMatrixMarketSystem(Teuchos::RCP<Teuchos::ParameterList>, 
	    Teuchos::RCP<CRS_MATRIX>             &A,
        Teuchos::RCP<MV>                     &b );
        
	Teuchos::RCP<CRS_MATRIX> applyShift(
    	Teuchos::RCP<CRS_MATRIX>,
    	Teuchos::RCP<Teuchos::ParameterList>);
     
    Teuchos::RCP<CRS_MATRIX> computeBlockDiagPrec(unsigned int);
                 
};

}

#endif // mc_solvers_LinearSystemFactory_hh

