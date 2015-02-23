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
#include "AdjointMcKernel.hh"
#include "ForwardMcKernel.hh"
#include "PolynomialFactory.hh"
#include "Kokkos_ExecPolicy.hpp"

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
 *  - weight_cutoff(SCALAR)   : >0.0 (1.0e-6)
 *  - verbosity(string)       : "none", ("low"), "medium", "high"
 */
//---------------------------------------------------------------------------//
MonteCarloSolver::MonteCarloSolver(Teuchos::RCP<const MATRIX> A,
                                   Teuchos::RCP<Teuchos::ParameterList> pl )
  : AleaSolver(A,pl)
{
    // Get Monte Carlo sublist
    Teuchos::RCP<Teuchos::ParameterList> mc_pl =
        Teuchos::sublist(pl,"Monte Carlo");

    // Override verbosity if present on sublist
    AleaSolver::setParameters(mc_pl);

    // Determine forward or adjoint
    std::string type = mc_pl->get("mc_type","adjoint");
    if( type == "forward" )
        d_type = FORWARD;
    else
        d_type = ADJOINT;

    // Get parameters off of PL
    std::string estimator = mc_pl->get<std::string>("estimator","expected_value");
    TEUCHOS_TEST_FOR_EXCEPT( estimator != "collision" &&
                             estimator != "expected_value" );
    d_use_expected_value = (estimator == "expected_value");

    if( d_type == FORWARD )
        d_use_expected_value = false;

    // Get number of requested threads
    d_num_threads = mc_pl->get<int>("num_threads",1);

    d_num_histories      = mc_pl->get<int>("num_histories",1000);
    d_weight_cutoff      = mc_pl->get<SCALAR>("weight_cutoff",1.0e-6);
    d_start_wt_factor    = mc_pl->get<SCALAR>("start_weight_factor",1.0);
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
    TEUCHOS_ASSERT( b_A != Teuchos::null );

    // Create Polynomial
    Teuchos::RCP<Polynomial> poly = PolynomialFactory::buildPolynomial(b_A,b_pl);
    TEUCHOS_ASSERT( poly != Teuchos::null );

    // Determine basis
    Teuchos::RCP<Teuchos::ParameterList> poly_pl =
        Teuchos::sublist(b_pl,"Polynomial");
    std::string basis_type = poly_pl->get("polynomial_basis","neumann");
    Teuchos::RCP<PolynomialBasis> basis( new PolynomialBasis(basis_type) );
    if( basis_type == "arbitrary" )
    {
        TEUCHOS_ASSERT( poly_pl->isType<SCALAR>("polynomial_basis_alpha") );
        TEUCHOS_ASSERT( poly_pl->isType<SCALAR>("polynomial_basis_beta") );
        SCALAR alpha = poly_pl->get<SCALAR>("polynomial_basis_alpha");
        SCALAR beta  = poly_pl->get<SCALAR>("polynomial_basis_beta");
        basis->setBasisCoefficients(alpha,beta);
    }

    // Get coefficients of polynomial in desired basis
    Teuchos::ArrayRCP<const SCALAR> coeffs = poly->getCoeffs(*basis);
    TEUCHOS_ASSERT( !coeffs.is_null() );
    Kokkos::resize(d_coeffs,coeffs.size());
    view_type::HostMirror coeffs_host = Kokkos::create_mirror_view(d_coeffs);
    std::copy(coeffs.begin(),coeffs.end(),&coeffs_host(0));
    Kokkos::deep_copy(d_coeffs,coeffs_host);

    // Get Monte Carlo sublist
    Teuchos::RCP<Teuchos::ParameterList> mc_pl =
        Teuchos::sublist(b_pl,"Monte Carlo");

    // Create Monte Carlo data
    d_mc_data = Teuchos::rcp(
        new MC_Data(b_A,basis,b_pl) );
    convertMatrices(d_mc_data->getIterationMatrix(),
                    d_mc_data->getProbabilityMatrix(),
                    d_mc_data->getWeightMatrix());

    if( d_num_histories%d_num_threads!= 0  && b_verbosity>=LOW )
    {
        std::cout << "WARNING: Requested number of histories ("
            << d_num_histories << ") is not divisible by the number "
            << "of threads (" << d_num_threads << "), ";
        d_num_histories = (d_num_histories/d_num_threads+1)*d_num_threads;
        std::cout << d_num_histories << " histories will be performed."
            << std::endl;
    }

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
    TEUCHOS_ASSERT( d_initialized );

    d_apply_count++;

    // For now we only support operating on a single vector
    TEUCHOS_ASSERT( x.getNumVectors() == 1 );

    // Check for early exit if x==0
    Teuchos::ArrayRCP<double> xnorm(1);
    x.norm2(xnorm());
    if( xnorm[0] == 0.0 )
    {
        y.putScalar(0.0);
        return;
    }

    LO N = x.getLocalLength();

    // Cache data components of x and y, on-the-fly access is SLOW
    const Teuchos::ArrayRCP<const SCALAR> x_data = x.getData(0);
    const Teuchos::ArrayRCP<SCALAR> y_data = y.getDataNonConst(0);

    int league_size = d_num_threads;
    int team_size = 1;
    Kokkos::TeamPolicy<DEVICE> exec_policy(league_size,team_size);

    if( d_type == FORWARD )
    {
        /*
        int histories_per_state = d_num_histories / N;

        if( d_num_histories % N != 0 )
            histories_per_state++;
        TEUCHOS_ASSERT( histories_per_state > 0 );

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
    else if( d_type == ADJOINT )
    {
        // Build initial probability and weight distributions
        // Should probably switch to doing this in a Kernel
        view_type start_cdf("starting_cdf",N);
        view_type start_wt("starting_wt",N);
        view_type::HostMirror start_cdf_host = Kokkos::create_mirror_view(start_cdf);
        view_type::HostMirror start_wt_host  = Kokkos::create_mirror_view(start_wt);
        for( LO i=0; i<N; ++i )
        {
            start_cdf_host(i) =
                SCALAR_TRAITS::pow(SCALAR_TRAITS::magnitude(x_data[i]),
                                   d_start_wt_factor);
        }
        SCALAR pdf_sum = std::accumulate(&start_cdf_host(0),&start_cdf_host(N-1)+1,0.0);
        TEUCHOS_ASSERT( pdf_sum > 0.0 );
        std::transform(&start_cdf_host(0),&start_cdf_host(N-1)+1,&start_cdf_host(0),
                       [pdf_sum](SCALAR x){return x/pdf_sum;});
        std::transform(x_data.begin(),x_data.end(),&start_cdf_host(0),
                       &start_wt_host(0),
                       [](SCALAR x, SCALAR y){return y==0.0 ? 0.0 : x/y;});
        std::partial_sum(&start_cdf_host(0),&start_cdf_host(N-1)+1,&start_cdf_host(0));
        Kokkos::deep_copy(start_cdf,start_cdf_host);
        Kokkos::deep_copy(start_wt,start_wt_host);

        int histories_per_thread = d_num_histories / d_num_threads;

        // Create temporary storage on device and mirror it on the host
        const view_type y_device("result",N);
        const view_type::HostMirror y_mirror =
            Kokkos::create_mirror_view(y_device);

        // Create kernel for performing group of MC histories
        std::cout << "Building AdjointMcKernel" << std::endl;
        AdjointMcKernel kernel(d_H,d_P,d_W,d_inds,d_offsets,d_coeffs,
                               start_cdf,start_wt,histories_per_thread,
                               d_use_expected_value,b_verbosity>=HIGH);

        // Execute Kokkos kernel on device
        std::cout << "Executing parallel_reduce" << std::endl;
        Kokkos::parallel_reduce( exec_policy, kernel, y_mirror);

        SCALAR scale_factor = 1.0 / static_cast<SCALAR>(d_num_histories);

        Kokkos::deep_copy(y_mirror,y_device);

        for( LO i=0; i<N; ++i )
        {
            y_data[i] = scale_factor*y_mirror(i);
        }

        // For expected value estimator, need to add C_0*x
        if( d_use_expected_value )
        {
            y.update(d_coeffs[0],x,1.0);
        }
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

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Convert MC matrices from CrsMatrix format to lightweight views.
 *
 * Data access in the Tpetra CrsMatrix implementation is not thread safe so
 * we get threadsafe views here to hand to the MC kernel.
 */
//---------------------------------------------------------------------------//
void MonteCarloSolver::convertMatrices(Teuchos::RCP<const MATRIX> H,
                                       Teuchos::RCP<const MATRIX> P,
                                       Teuchos::RCP<const MATRIX> W)
{
    TEUCHOS_ASSERT( H->isLocallyIndexed() );
    TEUCHOS_ASSERT( P->isLocallyIndexed() );
    TEUCHOS_ASSERT( W->isLocallyIndexed() );
    TEUCHOS_ASSERT( H->supportsRowViews() );
    TEUCHOS_ASSERT( P->supportsRowViews() );
    TEUCHOS_ASSERT( W->supportsRowViews() );

    LO numRows     = H->getNodeNumRows();
    GO numNonZeros = H->getNodeNumEntries();

    Kokkos::resize(d_H,numNonZeros);
    Kokkos::resize(d_P,numNonZeros);
    Kokkos::resize(d_W,numNonZeros);
    Kokkos::resize(d_inds,numNonZeros);
    Kokkos::resize(d_offsets,numRows+1);

    // Mirror views on host
    view_type::HostMirror H_host      = Kokkos::create_mirror_view(d_H);
    view_type::HostMirror P_host      = Kokkos::create_mirror_view(d_P);
    view_type::HostMirror W_host      = Kokkos::create_mirror_view(d_W);
    ord_view::HostMirror inds_host    = Kokkos::create_mirror_view(d_inds);
    ord_view::HostMirror offsets_host = Kokkos::create_mirror_view(d_offsets);

    Teuchos::ArrayView<const double> H_row, P_row, W_row;
    Teuchos::ArrayView<const int> ind_row;

    // Copy data from CrsMatrix into Kokkos host mirrors
    LO count = 0;
    for( LO irow=0; irow<numRows; ++irow )
    {
        // Get row views
        H->getLocalRowView(irow,ind_row,H_row);
        P->getLocalRowView(irow,ind_row,P_row);
        W->getLocalRowView(irow,ind_row,W_row);

        // Copy into host mirror
        std::copy( H_row.begin(), H_row.end(), &H_host(count) );
        std::copy( P_row.begin(), P_row.end(), &P_host(count) );
        std::copy( W_row.begin(), W_row.end(), &W_host(count) );
        std::copy( ind_row.begin(), ind_row.end(), &inds_host(count) );
        count += ind_row.size();
        offsets_host(irow+1) = count;
    }

    // Copy to device
    Kokkos::deep_copy(d_H,H_host);
    Kokkos::deep_copy(d_P,P_host);
    Kokkos::deep_copy(d_W,W_host);
    Kokkos::deep_copy(d_inds,inds_host);
    Kokkos::deep_copy(d_offsets,offsets_host);
}

} // namespace alea

