//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SyntheticAcceleration.hh
 * \author Steven Hamilton
 * \brief  Perform Richardson iteration.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_MultiSplitting_hh
#define mc_solvers_MultiSplitting_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "MatrixMarket_Tpetra.hpp"

#include "AleaTypedefs.hh"
#include "AleaSolver.hh"

#include "LinearSystem_MultiSplitting.hh"
#include "LinearSystemFactory.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class MultiSplitting
 * \brief Solver linear system using Multi-Splitting techniques
 */
//---------------------------------------------------------------------------//


class MultiSplitting 
{
  public:

    //! Enumeration describing level of output produced by solver
    enum Verbosity_Level
    {
        //! Produce no output.
        NONE,
        //! Minimal output, usually only a final summary.
        LOW,
        //! More output, commonly updates on the iteration process.
        MEDIUM,
        /*! Significant output, possibly extra information about each
            iteration. */
        HIGH,
        //! Highest level of output, may only be meaningful to developers.
        DEBUG
    };


    MultiSplitting(Teuchos::RCP<Teuchos::ParameterList> &);
         
    Teuchos::RCP<const CRS_MATRIX> getMatrix() const { return d_A; }     

    // Implementation of solver
    void solve(Teuchos::RCP<MV> &x);

  private:

    Teuchos::RCP<const CRS_MATRIX> d_A;
    Teuchos::RCP<const MV> d_b;

    // Parameter list
    Teuchos::RCP<Teuchos::ParameterList> d_pl;

    // Divergence tolerance
    MAGNITUDE d_divergence_tol;
    
    // Number of iterations employed
    mutable LO b_num_iters;
    
    // Maximal number of iterations allowed
    LO b_max_iterations;
 
    // Tolerance for the solver
    MAGNITUDE b_tolerance;
    
    // Verbosity level
    Verbosity_Level b_verbosity;

    // Inner solver type
    std::string d_inner_solver;

    // Numver of partitionings
    unsigned int d_num_blocks;
    
    // Parameters for the generations of subdomains
    Teuchos::RCP<LinearSystem_MultiSplitting> d_multisplitting;

    //Method to compute a MultiSplittng step
    Teuchos::RCP<MV> computeIteration();
    
};

}

#endif // mc_solvers_SyntheticAcceleration_hh

