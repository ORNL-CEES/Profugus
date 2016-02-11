//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MonteCarloEigenSolver.cc
 * \author Massimiliano Lupo Pasini
 * \brief  Perform Monte Carlo on eigenvalue problems
 */
//---------------------------------------------------------------------------//

#include <Alea/config.h>

#include <iterator>
#include <string>

#include "MonteCarloEigenSolver.hh"
#include "EigenMcAdaptive.hh"
//#include "ForwardMcKernel.hh"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "harness/DBC.hh"

//#ifdef USE_CUDA
//#include "EigenMcCuda.hh"
//#endif

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziMultiVecTraits.hpp"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param A  Problem matrix
 * \param pl ParameterList
 *
 * Behavior is controlled by following PL entries on the nested "Monte Carlo"
 * sublist:
 *  - mc_type(string)         : "forward" or ("adjoint")
 *  - estimator(string)       : "collision" or ("expected_value")
 *  - num_histories(int)      : >0 (1000)
 *  - verbosity(string)       : "none", ("low"), "medium", "high"
 */
//---------------------------------------------------------------------------//
MonteCarloEigenSolver::MonteCarloEigenSolver(Teuchos::RCP<const MATRIX> A,
                                   Teuchos::RCP<Teuchos::ParameterList> pl )
  : AleaSolver(A,pl)
  , d_rand_pool(pl->get("random_seed",31891))
{
    // Get Monte Carlo sublist
    d_mc_pl = Teuchos::sublist(pl,"Monte Carlo");

    // Override verbosity if present on sublist
    AleaSolver::setParameters(d_mc_pl);

    // Type of Kokkos kernel
    std::string kernel_type = d_mc_pl->get("kernel_type","adaptive");
    VALIDATE(kernel_type == "adaptive"        ||
             kernel_type == "cuda",
             "Invalid kernel_type.");

    if( kernel_type == "adaptive" )
        d_kernel_type = ADAPTIVE;
    else if( kernel_type == "cuda" )
        d_kernel_type = CUDA;

    d_num_histories = d_mc_pl->get<int>("num_histories",1000);
    d_init_count = 0;

    initialize();
    d_initialized = true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build data necessary for MC solve
 *
 * This function builds polynomial and MC data based on currently defined
 * matrices.  This is separate from the constructor to allow this object to
 * operate on a different matrix than the one that was used at construction.
 */
//---------------------------------------------------------------------------//
void MonteCarloEigenSolver::initialize()
{
    REQUIRE( b_A != Teuchos::null );

    // Create Monte Carlo data
    d_mc_data = Teuchos::rcp(new EigenMC_Data(b_A,b_pl));
    d_mc_data_kokkos = d_mc_data->createKokkosViews();

    b_label = "MonteCarloEigenSolver";
    d_initialized = true;
    d_init_count++;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform adjoint Monte Carlo process
 */
//---------------------------------------------------------------------------//
void MonteCarloEigenSolver::applyImpl(const MV &x, MV &y) const
{
    REQUIRE( d_initialized );

    d_apply_count++;

    // For now we only support operating on a single vector
    REQUIRE( x.getNumVectors() == 1 );

    // Check for early exit if x==0
    Teuchos::ArrayRCP<double> xnorm(1);
    x.norm2(xnorm());
    if( xnorm[0] == 0.0 )
    {
        y.putScalar(0.0);
        return;
    }

    LO N = x.getLocalLength();
    if( d_kernel_type == CUDA )
    {
/*#ifdef USE_CUDA
        EigenMcCuda solver(d_mc_data,d_mc_pl);

        solver.solve(x,y);
#else
        VALIDATE(false,"Cuda must be enabled to use Cuda kernel.");
#endif*/
    }    
    else if( d_kernel_type == ADAPTIVE )
    {
    	EigenMcAdaptive solver(d_mc_data,d_mc_pl,d_rand_pool);
    	solver.solve(x,y);
    }

    // There isn't a proper iteration count for MC
    // We use the number of histories as a stand-in
    b_num_iters = d_num_histories;
}

} // namespace alea

