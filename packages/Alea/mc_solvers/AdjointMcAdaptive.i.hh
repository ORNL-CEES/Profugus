//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdjointMcAdaptive.cc
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_AdjointMcAdaptive_i_hh
#define mc_solver_AdjointMcAdaptive_i_hh

#include <iterator>
#include <random>
#include <cmath>

#include "AdjointMcAdaptive.hh"
#include "MC_Components.hh"
#include "utils/String_Functions.hh"
#include "harness/Warnings.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param P Views into entries of probability matrix
 * \param W Views into entries of weight matrix
 * \param inds Views into nonzeros indices
 * \param offsets Starting indices for each matrix row
 * \param coeffs Polynomial coefficients
 * \param pl Problem parameters
 */
//---------------------------------------------------------------------------//
AdjointMcAdaptive::AdjointMcAdaptive(Teuchos::RCP<const MC_Data> mc_data,
        Teuchos::RCP<Teuchos::ParameterList> pl,
        generator_pool rand_pool)

  : d_N(mc_data->getIterationMatrix()->getGlobalNumRows())
  , d_rand_pool(rand_pool)
  , d_rand_gen(d_rand_pool.get_state())
{
    // Get parameters
    d_num_histories      = pl->get("num_histories",1000);
    d_max_history_length = pl->get("max_history_length",1000);
    d_weight_cutoff      = pl->get("weight_cutoff",0.0);

    // Determine type of tally
    std::string estimator = pl->get<std::string>("estimator",
                                                 "expected_value");
    VALIDATE(estimator == "collision" ||
             estimator == "expected_value",
             "Only collision and expected_value estimators are available.");
    d_use_expected_value = (estimator == "expected_value");

    // Should we print anything to screen
    d_verbosity = profugus::to_lower(pl->get("verbosity","low"));

    extractMatrices(mc_data);
}

//---------------------------------------------------------------------------//
// Solve problem using Monte Carlo
//---------------------------------------------------------------------------//
void AdjointMcAdaptive::solve(const MV &b, MV &x)
{
    // Containers to hold starting cdf and wt arrays
    Teuchos::ArrayRCP<double> start_cdf, start_wt;

    // Build initial probability and weight distributions
    Teuchos::ArrayRCP<const double> b_data = b.getData(0);
    build_initial_distribution(b_data,start_cdf,start_wt);

    Teuchos::ArrayRCP<double> x_data = x.getDataNonConst(0);

    // Storage for current row of H, P, W
    Teuchos::ArrayView<const double> h_row, p_row, w_row;
    Teuchos::ArrayView<const int> ind_row;

    int state = -1;
    double wt = 0.0;
    double init_wt = 0.0;
    for( int i=0; i<d_num_histories; ++i )
    {
        // Get initial state for this history by sampling from start_cdf
        initializeHistory(state,wt,start_cdf,start_wt,h_row,p_row,w_row,ind_row);
        init_wt = wt;

        // With expected value estimator we start on stage 1 because
        // zeroth order term is added explicitly at the end
        int stage = d_use_expected_value ? 1 : 0;

        // Perform initial tally
        tallyContribution(state,wt,x_data,h_row,ind_row);

        for( ; stage<=d_max_history_length; ++stage )
        {
            // Move to new state
            getNewState(state,wt,h_row,p_row,w_row,ind_row);

            // Tally
            tallyContribution(state,wt,x_data,h_row,ind_row);

            // Check weight cutoff
            if( std::abs(wt/init_wt) < d_weight_cutoff )
                break;
        }
    }

    // Normalize by number of histories
    double scale_factor = 1.0 / static_cast<double>(d_num_histories);
    std::transform(x_data.begin(),x_data.end(),x_data.begin(),
        [scale_factor](double val){return val*scale_factor;});

    // Add rhs for expected value
    if( d_use_expected_value )
        x.update(1.0,b,1.0);
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Extract matrices into ArrayView objects for faster data access
//---------------------------------------------------------------------------//
void AdjointMcAdaptive::extractMatrices(Teuchos::RCP<const MC_Data> mc_data)
{
    d_H.resize(d_N);
    d_P.resize(d_N);
    d_W.resize(d_N);
    d_ind.resize(d_N);

    Teuchos::RCP<const MATRIX> H = mc_data->getIterationMatrix();
    Teuchos::RCP<const MATRIX> P = mc_data->getProbabilityMatrix();
    Teuchos::RCP<const MATRIX> W = mc_data->getWeightMatrix();

    Teuchos::ArrayView<const double> val_row;
    Teuchos::ArrayView<const int>    ind_row;
    for( int i=0; i<d_N; ++i )
    {
        // Extract row i of matrix
        H->getLocalRowView(i,ind_row,val_row);
        d_H[i] = val_row;
        P->getLocalRowView(i,ind_row,val_row);
        d_P[i] = val_row;
        W->getLocalRowView(i,ind_row,val_row);
        d_W[i] = val_row;
        d_ind[i] = ind_row;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//
void AdjointMcAdaptive::tallyContribution(int state, double wt,
        Teuchos::ArrayRCP<double>        x,
        Teuchos::ArrayView<const double> h_row,
        Teuchos::ArrayView<const int>    ind_row ) const
{
    if( d_use_expected_value )
    {
        REQUIRE( h_row.size() == ind_row.size() );

        // For expected value estimator, loop over current row and add
        // contributions corresponding to each element
        for( int i=0; i<h_row.size(); ++i )
        {
            x[ind_row[i]] += wt * h_row[i];
        }
    }
    else
    {
        // Collision estimator just adds weight
        x[state] += wt;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize history into new state
 */
//---------------------------------------------------------------------------//
void AdjointMcAdaptive::initializeHistory(int &state, double &wt,
        Teuchos::ArrayRCP<const double>   start_cdf,
        Teuchos::ArrayRCP<const double>   start_wt,
        Teuchos::ArrayView<const double> &h_row,
        Teuchos::ArrayView<const double> &p_row,
        Teuchos::ArrayView<const double> &w_row,
        Teuchos::ArrayView<const int>    &ind_row)
{
    // Generate random number
    double rand = Kokkos::rand<generator_type,double>::draw(d_rand_gen);

    // Sample cdf to get new state
    auto elem = std::lower_bound(start_cdf.begin(),start_cdf.end(),rand);

    if( elem == start_cdf.end() )
    {
        // Invalidate all row data
        h_row   = Teuchos::ArrayView<const double>();
        p_row   = Teuchos::ArrayView<const double>();
        w_row   = Teuchos::ArrayView<const double>();
        ind_row = Teuchos::ArrayView<const int>();
        return;
    }

    // Modify weight and update state
    state = elem-start_cdf.begin();
    wt    = start_wt[state];

    // Get new rows for this state
    h_row   = d_H[state];
    p_row   = d_P[state];
    w_row   = d_W[state];
    ind_row = d_ind[state];
}
//---------------------------------------------------------------------------//
/*!
 * \brief Get new state by sampling from cdf
 */
//---------------------------------------------------------------------------//
void AdjointMcAdaptive::getNewState(int &state, double &wt,
        Teuchos::ArrayView<const double> &h_row,
        Teuchos::ArrayView<const double> &p_row,
        Teuchos::ArrayView<const double> &w_row,
        Teuchos::ArrayView<const int>    &ind_row)
{
    // Generate random number
    double rand = Kokkos::rand<generator_type,double>::draw(d_rand_gen);

    // Sample cdf to get new state
    auto elem = std::lower_bound(p_row.begin(),p_row.end(),rand);

    if( elem == p_row.end() )
    {
        // Invalidate all row data
        h_row   = Teuchos::ArrayView<const double>();
        p_row   = Teuchos::ArrayView<const double>();
        w_row   = Teuchos::ArrayView<const double>();
        ind_row = Teuchos::ArrayView<const int>();
        return;
    }

    // Modify weight and update state
    int index = elem - p_row.begin();
    state  =  ind_row[index];
    wt    *=  w_row[index];

    // Get new rows for this state
    h_row   = d_H[state];
    p_row   = d_P[state];
    w_row   = d_W[state];
    ind_row = d_ind[state];
}

//---------------------------------------------------------------------------//
// Build initial cdf and weights
//---------------------------------------------------------------------------//
void AdjointMcAdaptive::build_initial_distribution(
        Teuchos::ArrayRCP<const double> b,
        Teuchos::ArrayRCP<double>      &cdf,
        Teuchos::ArrayRCP<double>      &wt) const
{
    cdf.resize(d_N);
    wt.resize(d_N);

    // First take absolute value of b
    for( int i=0; i<d_N; ++i )
    {
        cdf[i] = std::abs(b[i]);
    }

    // Normalize to get a PDF
    double pdf_sum = std::accumulate(cdf.begin(),cdf.end(),0.0);
    ENSURE( pdf_sum > 0.0 );
    std::transform(cdf.begin(),cdf.end(),cdf.begin(),
            [pdf_sum](double val){return val/pdf_sum;});

    // Compute weight vector s.t. b = p * wt
    std::transform(b.begin(),b.end(),cdf.begin(),wt.begin(),
            [](double u, double v){return v==0.0 ? 0.0 : u/v;});

    // Convert PDF to CDF
    std::partial_sum(cdf.begin(),cdf.end(),cdf.begin());
}

} // namespace alea

#endif // mc_solver_AdjointMcAdaptive_i_hh
