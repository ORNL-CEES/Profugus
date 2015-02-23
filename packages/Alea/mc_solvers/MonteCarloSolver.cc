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
    d_mc_pl = Teuchos::sublist(pl,"Monte Carlo");

    // Override verbosity if present on sublist
    AleaSolver::setParameters(d_mc_pl);

    // Determine forward or adjoint
    std::string type = d_mc_pl->get("mc_type","adjoint");
    if( type == "forward" )
        d_type = FORWARD;
    else
        d_type = ADJOINT;

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

    // Create Monte Carlo data
    d_mc_data = Teuchos::rcp(
        new MC_Data(b_A,basis,b_pl) );
    convertMatrices(d_mc_data->getIterationMatrix(),
                    d_mc_data->getProbabilityMatrix(),
                    d_mc_data->getWeightMatrix());

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
        // Create kernel for performing group of MC histories
        AdjointMcKernel kernel(d_H,d_P,d_W,d_inds,d_offsets,d_coeffs,d_mc_pl);

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

