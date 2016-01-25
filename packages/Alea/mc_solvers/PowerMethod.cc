//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PowerMethod.cc
 * \author Massimiliano Lupo Pasni
 * \brief  Perform power iteration on eigenvalue system
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "PowerMethod.hh"
#include "EigenSolverFactory.hh"
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
 * Behavior is controlled by following PL entries on the nested "power iteration"
 * sublist:
 *  - max_iterations(int)    : >0 (1000)
 *  - tolerance(SCALAR)      : >0.0 (1.0e-6)
 *  - verbosity(string)      : "none", ("low"), "medium", "high"
 */
//---------------------------------------------------------------------------//
PowerMethod::PowerMethod(
        Teuchos::RCP<const MATRIX> A,
        Teuchos::RCP<Teuchos::ParameterList> pl )
  : AleaSolver(A,pl)
{
    // Get belos pl
    Teuchos::RCP<Teuchos::ParameterList> power_pl =
        Teuchos::sublist(pl,"PowerMethod");

    d_divergence_tol = power_pl->get("divergence_tolerance",1.0e4);

    // Override default parameters if present on sublist
    this->setParameters(power_pl);

    b_label = "PowerMethod";
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Perform Richardson iteration.
 */
//---------------------------------------------------------------------------//
void PowerMethod::applyImpl(const MV &x, MV &y) const
{
    // For now we only support operating on a single vector
    REQUIRE( x.getNumVectors() == 1 );

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1, 2);

    SCALAR eig_old = dis(gen);
    SCALAR eig_new = 0.0;

    // Compute initial residual
    MV r(y.getMap(),1);
    r.update(1.0,x,0.0);
    y.update(1.0,x,0.0);
    b_num_iters = 0;

    Teuchos::ArrayRCP<SCALAR> r_data = r.getDataNonConst(0);
    Teuchos::ArrayRCP<SCALAR> y_data = y.getDataNonConst(0);

    Teuchos::ArrayRCP<MAGNITUDE> y_norm(1);
    y.norm2(y_norm());
    unsigned int N = y_data.size();
 
    while( true )
    {
        b_num_iters++;

        // Check convergence on true (rather than preconditioned) residual
        y.norm2(y_norm());

    	for(unsigned int i = 0; i!=N; ++i)
		y_data[i] /= y_norm[0];

        // Compute residual r = A*y
        b_A->apply(y,r);

        eig_new = 0.0;

        for(unsigned int i = 0; i!=N; ++i)
		eig_new += r_data[i]*y_data[i];

	std::cout<<eig_new<<std::endl;

        if( b_verbosity >= HIGH )
        {
            std::cout << "Approximated eigenvalue at iteration " << b_num_iters
                << " is " << eig_new << std::endl;
        }

        // Check for convergence
        if( std::abs(eig_old - eig_new)/std::abs(eig_old) < b_tolerance )
        {
            if( b_verbosity >= LOW )
            {
                std::cout << "Power Iteration converged after "
                    << b_num_iters << " iterations." << std::endl;
            }
            break;
        }

        // Check for max iterations
        if( b_num_iters >= b_max_iterations )
        {
            if( b_verbosity >= LOW )
            {
                std::cout << "Power Iteration reached maximum iteration "
                    << "count with relative error of "
                    << std::abs(eig_old - eig_new)/std::abs(eig_old) << std::endl;
            }
            b_num_iters = -1;
            break;
        }

        // Check for divergence
        if( std::abs(eig_old - eig_new)/std::abs(eig_old) > d_divergence_tol)
        {
            if( b_verbosity >= LOW )
            {
                std::cout << "Richardson Iteration diverged after "
                    << b_num_iters << " iterations." << std::endl;
            }
            b_num_iters = -b_num_iters;
            break;
        }
        
        eig_old = eig_new;
        y.update(1.0,r,0.0);
    }

}

} // namespace alea

