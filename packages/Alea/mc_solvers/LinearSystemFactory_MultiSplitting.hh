//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSystemFactory_MultiSplitting.hh
 * \author Steven Hamilton
 * \brief  Construct EpetraCrsMatrix from ParameterList.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_LinearSystemFactory_MultiSplitting_hh
#define mc_solvers_LinearSystemFactory_MultiSplitting_hh

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
class LinearSystemFactory_MultiSplitting
{
  private:

    // Pure static, disallow construction
    LinearSystemFactory_MultiSplitting(){};

  public:

    static void LinearSystemFactory_MultiSplitting::createPartitions( 
          Teuchos::RCP<Teuchos::ParameterList>);    

    static splitting buildSplitting(unsigned int);
        
    static void LinearSystemFactory_MultiSplitting::buildSystem( 
        Teuchos::RCP<Teuchos::ParameterList>);    

  private:

    Teuchos::RCP<CRS_MATRIX> d_A;
    Teuchos::RCP<MV> d_b; 
    
    Teuchos::ArrayRCP< ENDPOINTS > d_partitions;

	static void LinearSystemFactory_MultiSplitting::buildMatrixMarketSystem(
        Teuchos::RCP<Teuchos::ParameterList>,
        Teuchos::RCP<CRS_MATRIX> &,
        Teuchos::RCP<MV>         &);
        
	static Teuchos::RCP<CRS_MATRIX> LinearSystemFactory_MultiSplitting::applyShift(
    	Teuchos::RCP<CRS_MATRIX>,
    	Teuchos::RCP<Teuchos::ParameterList>);
     

    static Teuchos::RCP<CRS_MATRIX> 
        LinearSystemFactory_MultiSplitting::computeBlockDiagPrec(unsigned int);
                
    unsigned int d_num_blocks;
    SCALAR d_overlap;
    std::string d_inner_solver;    
};

}

#endif // mc_solvers_LinearSystemFactory_hh

