//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ForwardMcAdaptive.cc
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_ForwardMcAdaptive_i_hh
#define mc_solver_ForwardMcAdaptive_i_hh

#include <iterator>
#include <random>
#include <cmath>

#include "ForwardMcAdaptive.hh"
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
ForwardMcAdaptive::ForwardMcAdaptive(
        Teuchos::RCP<const MC_Data> mc_data,
        const const_scalar_view     coeffs,
        Teuchos::RCP<Teuchos::ParameterList> pl,
        generator_pool rand_pool)

  : d_N(mc_data->getIterationMatrix()->getGlobalNumRows())
  , d_coeffs(coeffs) // Modified by Max
  , d_rand_pool(rand_pool)
  , d_rand_gen(d_rand_pool.get_state())
{
    // Get parameters
    d_max_num_histories  = pl->get("num_histories",1000);
    d_max_history_length = d_coeffs.size()-1;
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
void ForwardMcAdaptive::solve(const MV &b, MV &x)
{
    std::cout<<"Forward solver called"<<std::endl;
    Teuchos::ArrayRCP<double> x_data = x.getDataNonConst(0);
    std::fill( x_data.begin(), x_data.end(), 0.0 );

    // Storage for current row of H, P, W
    Teuchos::ArrayView<const double> h_row, p_row, w_row;
    Teuchos::ArrayView<const int> ind_row;

    double x_batch;
    Teuchos::ArrayRCP<double> variance(d_N);
    double variance_batch;

    Teuchos::ArrayRCP<const double> b_data = b.getData(0);

    int state = -1;
    double wt = 0.0;

    int total_histories = 0;

    for (int entry=0; entry < d_N; ++entry)
    {
    	double rel_std_dev = 1e6;
    	int batch=0;
    	int num_histories = 0;
    	batch=0;

        while( ( rel_std_dev > d_tolerance && num_histories < d_max_num_histories ) ||
         ( rel_std_dev == 0.0 && num_histories < d_max_num_histories ) )
        {
            batch++;
            x_batch = 0.0;
            variance_batch = 0.0;

            for( int i=0; i<d_batch_size; ++i )
            {
                double x_history = 0.0;
                int stage = 0 ;
                int init_wt = 1.0;
                wt = 1.0 ;

                // Perform initial tally
                tallyContribution(d_coeffs[stage]*wt*b_data[entry],x_history);

                // Get new rows for this state
                h_row   = d_H[entry];
                p_row   = d_P[entry];
                w_row   = d_W[entry];
                ind_row = d_ind[entry];

                for( ; stage<=d_max_history_length; ++stage )
                {
                    // Move to new state
                    getNewState(state,wt,h_row,p_row,w_row,ind_row);
                    if( h_row.size() == 0 )
                        break;

                    // Tally
                    tallyContribution(d_coeffs[stage]*wt*b_data[state],x_history);

                    // Check weight cutoff
                    if( std::abs(wt) < d_weight_cutoff )
                        break;
                }

                x_batch  += x_history;
            }

            variance_batch += x_batch * x_batch;

            // From the old variance, compute the new second moment
            variance[entry] = (variance[entry] * static_cast<double>(num_histories-1) +
                    x_data[entry]*x_data[entry]*static_cast<double>(num_histories) +
                    variance_batch);

            // Compute new mean
            x_data[entry] = (x_data[entry] * static_cast<double>(num_histories) +
                    x_batch) / static_cast<double>(num_histories+d_batch_size);

            num_histories += d_batch_size;

            variance[entry] = (variance[entry] - x_data[entry]*x_data[entry]*
                    static_cast<double>(num_histories)) /
                static_cast<double>(num_histories-1);

            if(x_data[entry] == 0.0)
            {
                rel_std_dev=0.0;
                break;
            }

            // Compute 1-norm of solution and variance of mean
            double std_dev = 0;
            double var = variance[entry] / static_cast<double>(num_histories);
            if( var > 0.0 )
                std_dev += std::sqrt(var);

            CHECK( static_cast<double>(std::abs(x_data[entry])) > 0.0 );
            rel_std_dev = static_cast<double>( std_dev / static_cast<double>(std::abs(x_data[entry])) );

        }
        std::cout << "Entry " << entry << " performed " << num_histories << " histories" << " with final std dev of " << rel_std_dev << std::endl;
        total_histories += num_histories;
	//std::cout<<rel_std_dev<<std::endl;
    }

    std::cout << "Performed " << total_histories << " total histories, "
        << " average of " <<
        static_cast<double>(total_histories)/static_cast<double>(d_N)
        << " per entry" << std::endl;

	double sol_1norm = 0.0;
	for(int i=0; i<d_N; ++i)
		sol_1norm += static_cast<double>( std::abs(x_data[i]) );

	double std_dev_1norm = 0.0;
	for(int i=0; i<d_N; ++i)
		std_dev_1norm += static_cast<double>( std::sqrt( variance[i] ) );


	if( d_verbosity == HIGH )
	{
	    std::cout << "Relative std dev: "  << static_cast<double>( std_dev_1norm / std_dev_1norm) << std::endl;
	}


/*    if( d_verbosity >= LOW )
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
    }*/
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Extract matrices into ArrayView objects for faster data access
//---------------------------------------------------------------------------//
void ForwardMcAdaptive::extractMatrices(Teuchos::RCP<const MC_Data> mc_data)
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
void ForwardMcAdaptive::tallyContribution(double wt, double& x)
{
        x += wt;
}

/*!
 * \brief Get new state by sampling from cdf
 */
//---------------------------------------------------------------------------//
void ForwardMcAdaptive::getNewState(int &state, double &wt,
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


} // namespace alea

#endif // mc_solver_AdjointMcAdaptive_i_hh
