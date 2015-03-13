//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MonteCarloSolver.cc
 * \author Steven Hamilton
 * \brief  Perform Adjoint Monte Carlo on linear system
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "MonteCarloSolver.hh"
#include "AdjointMcParallelFor.hh"
#include "AdjointMcParallelReduce.hh"
#include "AdjointMcEventKernel.hh"
//#include "ForwardMcKernel.hh"
#include "PolynomialFactory.hh"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "harness/DBC.hh"

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
MonteCarloSolver::MonteCarloSolver(Teuchos::RCP<const MATRIX> A,
                                   Teuchos::RCP<Teuchos::ParameterList> pl )
  : AleaSolver(A,pl)
  , d_rand_pool(pl->get("random_seed",31891))
{
    // Get Monte Carlo sublist
    d_mc_pl = Teuchos::sublist(pl,"Monte Carlo");

    // Override verbosity if present on sublist
    AleaSolver::setParameters(d_mc_pl);

    // Determine forward or adjoint
    std::string mc_type = d_mc_pl->get("mc_type","adjoint");
    VALIDATE(mc_type == "adjoint" || mc_type == "forward",
             "Invalid mc_type.");
    if( mc_type == "forward" )
        d_mc_type = FORWARD;
    else
        d_mc_type = ADJOINT;

    // Type of Kokkos kernel
    std::string kernel_type = d_mc_pl->get("kernel_type","parallel_reduce");
    VALIDATE(kernel_type == "parallel_reduce" ||
             kernel_type == "parallel_for"    ||
             kernel_type == "event",
             "Invalid kernel_type.");
    if( kernel_type == "parallel_for" )
        d_kernel_type = PARALLEL_FOR;
    else if( kernel_type == "parallel_reduce" )
        d_kernel_type = PARALLEL_REDUCE;
    else if( kernel_type == "event" )
        d_kernel_type = EVENT;

    d_num_histories = d_mc_pl->get<int>("num_histories",1000);
    d_init_count = 0;
    d_initialized = false;
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
void MonteCarloSolver::initialize()
{
    REQUIRE( b_A != Teuchos::null );

    // Create Polynomial
    Teuchos::RCP<Polynomial> poly = PolynomialFactory::buildPolynomial(b_A,b_pl);
    REQUIRE( poly != Teuchos::null );

    // Determine basis
    Teuchos::RCP<Teuchos::ParameterList> poly_pl =
        Teuchos::sublist(b_pl,"Polynomial");
    std::string basis_type = poly_pl->get("polynomial_basis","neumann");
    Teuchos::RCP<PolynomialBasis> basis( new PolynomialBasis(basis_type) );
    if( basis_type == "arbitrary" )
    {
        CHECK( poly_pl->isType<SCALAR>("polynomial_basis_alpha") );
        CHECK( poly_pl->isType<SCALAR>("polynomial_basis_beta") );
        SCALAR alpha = poly_pl->get<SCALAR>("polynomial_basis_alpha");
        SCALAR beta  = poly_pl->get<SCALAR>("polynomial_basis_beta");
        basis->setBasisCoefficients(alpha,beta);
    }

    // Get coefficients of polynomial in desired basis
    Teuchos::ArrayRCP<const SCALAR> coeffs = poly->getCoeffs(*basis);
    CHECK( !coeffs.is_null() );
    Kokkos::resize(d_coeffs,coeffs.size());
    scalar_view::HostMirror coeffs_host = Kokkos::create_mirror_view(d_coeffs);
    std::copy(coeffs.begin(),coeffs.end(),&coeffs_host(0));
    Kokkos::deep_copy(d_coeffs,coeffs_host);

    // Create Monte Carlo data
    Teuchos::RCP<MC_Data> mc_data( new MC_Data(b_A,basis,b_pl) );
    d_mc_data = mc_data->createKokkosViews();

    b_label = "MonteCarloSolver";
    d_initialized = true;
    d_init_count++;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform adjoint Monte Carlo process
 */
//---------------------------------------------------------------------------//
void MonteCarloSolver::applyImpl(const MV &x, MV &y) const
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

    if( d_mc_type == FORWARD )
    {
        /*
        int histories_per_state = d_num_histories / N;

        if( d_num_histories % N != 0 )
            histories_per_state++;
        CHECK( histories_per_state > 0 );

        // Construct and execute Kokkos kernel
        ForwardMcKernel kernel(d_P(),d_W(),d_inds(),d_coeffs(),N,x_data(),
                               y_data(),histories_per_state,
                               b_verbosity>=HIGH);

        // Create kernel for performing group of MC histories
        Kokkos::parallel_for( exec_policy, kernel );

        SCALAR scale_factor = static_cast<SCALAR>(N) /
                              static_cast<SCALAR>(d_num_histories);
        std::transform(y_data.begin(),y_data.end(),y_data.begin(),
                       [scale_factor](SCALAR x){return x*scale_factor;});
                       */

    }
    else if( d_mc_type == ADJOINT && d_kernel_type == PARALLEL_REDUCE )
    {
        // Create kernel for performing group of MC histories
        AdjointMcParallelReduce kernel(d_mc_data,d_coeffs,d_mc_pl);

        kernel.solve(x,y);
    }
    else if( d_mc_type == ADJOINT && d_kernel_type == PARALLEL_FOR )
    {
        // Create kernel for performing group of MC histories
        AdjointMcParallelFor kernel(d_mc_data,d_coeffs,d_mc_pl);

        kernel.solve(x,y);
    }
    else if( d_mc_type == ADJOINT && d_kernel_type == EVENT )
    {
        // Create kernel for performing group of MC histories
        AdjointMcEventKernel kernel(d_mc_data,d_coeffs,d_mc_pl,d_rand_pool);

        kernel.solve(x,y);
    }

    if( b_verbosity >= LOW )
    {
        std::cout << "Performed " << d_num_histories
            << " histories" << std::endl;
    }

    // There isn't a proper iteration count for MC
    // We use the number of histories as a stand-in
    b_num_iters = d_num_histories;
}

} // namespace alea

