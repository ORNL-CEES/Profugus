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
AdjointMcAdaptive::AdjointMcAdaptive(
        Teuchos::RCP<const MC_Data> mc_data,
        Teuchos::RCP<Teuchos::ParameterList> pl,
        generator_pool rand_pool)

  : d_N(mc_data->getIterationMatrix()->getGlobalNumRows())
  , d_rand_pool(rand_pool)
  , d_rand_gen(d_rand_pool.get_state())
{
    // Get parameters
    d_max_num_histories  = pl->get("num_histories",1000);
    d_max_history_length = pl->get("max_history_length",1000);
    d_weight_cutoff      = pl->get("weight_cutoff",0.0);
    d_batch_size         = pl->get("batch_size",100);
    d_tolerance          = pl->get("tolerance",0.01);

    // Determine type of tally
    std::string estimator = pl->get<std::string>("estimator",
                                                 "expected_value");
    VALIDATE(estimator == "collision" ||
             estimator == "expected_value",
             "Only collision and expected_value estimators are available.");
    d_use_expected_value = (estimator == "expected_value");

    // Should we print anything to screen
    std::string verb = profugus::to_lower(pl->get("verbosity","low"));
    if( verb == "none" )
        d_verbosity = NONE;
    else if( verb == "low" )
        d_verbosity = LOW;
    else if( verb == "high" )
        d_verbosity = HIGH;

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
    std::fill( x_data.begin(), x_data.end(), 0.0 );

    // Storage for current row of H, P, W
    Teuchos::ArrayView<const double> h_row, p_row, w_row;
    Teuchos::ArrayView<const int> ind_row;

    Teuchos::ArrayRCP<double> x_history(d_N);
    Teuchos::ArrayRCP<double> x_batch(d_N);
    Teuchos::ArrayRCP<double> variance(d_N);
    Teuchos::ArrayRCP<double> variance_batch(d_N);

    int state = -1;
    double wt = 0.0;
    double init_wt = 0.0;

    double rel_std_dev = 1e6;
    int batch=0;
    int num_histories = 0;
    while( rel_std_dev > d_tolerance && num_histories < d_max_num_histories )
    {
        batch++;
        std::fill( x_batch.begin(), x_batch.end(), 0.0 );
        std::fill( variance_batch.begin(), variance_batch.end(), 0.0 );

        for( int i=0; i<d_batch_size; ++i )
        {
            std::fill( x_history.begin(), x_history.end(), 0.0 );

            // Get initial state for this history by sampling from start_cdf
            initializeHistory(state,wt,start_cdf,start_wt,h_row,p_row,w_row,
                ind_row);
            CHECK( h_row.size() > 0 );
            init_wt = wt;

            // With expected value estimator we start on stage 1 because
            // zeroth order term is added explicitly at the end
            int stage = d_use_expected_value ? 1 : 0;

            // Perform initial tally
            tallyContribution(state,wt,x_history,h_row,ind_row);

            for( ; stage<=d_max_history_length; ++stage )
            {
                // Move to new state
                getNewState(state,wt,h_row,p_row,w_row,ind_row);
                if( h_row.size() == 0 )
                    break;

                // Tally
                tallyContribution(state,wt,x_history,h_row,ind_row);

                // Check weight cutoff
                if( std::abs(wt/init_wt) < d_weight_cutoff )
                    break;
            }

            for( int i=0; i<d_N; ++i )
            {
                x_batch[i]  += x_history[i];
                variance_batch[i] += x_history[i]*x_history[i];
            }

        }

        // Compute new second moment and mean from batch results
        for( int i=0; i<d_N; ++i )
        {
            // From the old variance, compute the new second moment
            variance[i] = (variance[i] * static_cast<double>(num_histories-1) +
                x_data[i]*x_data[i]*static_cast<double>(num_histories) +
                variance_batch[i]);

            // Compute new mean
            x_data[i] = (x_data[i] * static_cast<double>(num_histories) +
                x_batch[i]) / static_cast<double>(num_histories+d_batch_size);
        }

        num_histories += d_batch_size;

        // Add rhs for expected value
        if( d_use_expected_value )
            x.update(1.0,b,1.0);

        // Subtract square of mean from second moment to get variance
        for( int i=0; i<d_N; ++i )
        {
            variance[i] = (variance[i] - x_data[i]*x_data[i]*
                static_cast<double>(num_histories)) /
                static_cast<double>(num_histories-1);
        }

        // Compute 1-norm of solution and variance of mean
        double soln_1norm = 0.0;
        double std_dev_1norm = 0.0;
        for( int i=0; i<d_N; ++i )
        {
            soln_1norm += std::abs(x_data[i]);
            double var = variance[i] / static_cast<double>(num_histories);
            if( var > 0.0 )
                std_dev_1norm += std::sqrt(var);
        }
        CHECK( soln_1norm > 0.0 );
        rel_std_dev = std_dev_1norm / soln_1norm;

        if( d_verbosity == HIGH )
        {
            std::cout << "Relative std dev after " << num_histories <<
                " histories : " << rel_std_dev << std::endl;
        }
    }

    if( d_verbosity >= LOW )
    {
        if( rel_std_dev < d_tolerance )
        {
            std::cout << "Converged with relative std dev of " << rel_std_dev
                << " using " << num_histories << " histories"
                << std::endl;
        }
        else
        {
            std::cout << "Did not converge, final relative std dev is " <<
                rel_std_dev << std::endl;
        }
    }
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
