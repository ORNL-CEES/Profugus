//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/ChebyshevIteration.cc
 * \author Steven Hamilton
 * \brief  Perform Adjoint Monte Carlo on linear system
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "ChebyshevIteration.hh"
#include "ChebyshevPolynomial.hh"
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
 * Behavior is controlled by following PL entries on the nested "Chebyshev"
 * sublist:
 *  max_iterations(int) : >0 (1000)
 *  tolerance(SCALAR)   : >0.0 (1.0e-6)
 *  verbosity(string)   : "none", ("low"), "medium", "high"
 */
//---------------------------------------------------------------------------//
ChebyshevIteration::ChebyshevIteration(
        Teuchos::RCP<const MATRIX> A,
        Teuchos::RCP<Teuchos::ParameterList> pl )
  : AleaSolver(A,pl)
{
    // Get belos pl
    Teuchos::RCP<Teuchos::ParameterList> cheb_pl =
        Teuchos::sublist(pl,"Chebyshev");

    // Override default parameters if present on sublist
    this->setParameters(cheb_pl);

    // Set eigenvalue parameters
    Teuchos::RCP<Teuchos::ParameterList> poly_pl =
        Teuchos::sublist(pl,"Polynomial");
    poly_pl->get("compute_eigenvalues",true); // "get" sets default value

    Teuchos::RCP<ChebyshevPolynomial> poly( new ChebyshevPolynomial(A,pl) );
    SCALAR lambda_min = poly->getLambdaMin();
    SCALAR lambda_max = poly->getLambdaMax();
    INSIST( !SCALAR_TRAITS::isnaninf(lambda_min),
            "Minimum eigenvalue is NaN or inf" );
    INSIST( !SCALAR_TRAITS::isnaninf(lambda_max),
            "Maximum eigenvalue is NaN or inf" );

    d_c = (lambda_max+lambda_min)/2.0;
    d_d  = (lambda_max-lambda_min)/2.0;

    b_label = "ChebyshevIteration";
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Perform Chebyshev iteration.
 */
//---------------------------------------------------------------------------//
void ChebyshevIteration::applyImpl(const MV &x, MV &y) const
{
    // For now we only support operating on a single vector
    REQUIRE( x.getNumVectors() == 1 );

    // Compute initial residual
    MV r(y.getMap(),1);
    r.update(1.0,x,0.0);

    Teuchos::ArrayRCP<MAGNITUDE> r_norm(1);
    r.norm2(r_norm());
    if( r_norm[0] == 0.0 )
    {
        if( b_verbosity >= LOW )
        {
            std::cout << "Input vector has zero norm."
                << "  Setting output vector to zero without iterating."
                << std::endl;
        }
        return;
    }

    MAGNITUDE r0_norm = r_norm[0];

    // Local vectors
    bool init_to_zero = true;
    MV delta(y.getMap(),y.getNumVectors(),init_to_zero);

    SCALAR a=0.0, b=-d_c;;

    b_num_iters = 0;
    while( true )
    {
        b_num_iters++;

        if( b_num_iters > 1 )
        {
            a = d_d*d_d / (4.0*b);
            b = -(d_c + a);
        }

        // delta = (1/b)*(-r + a*delta)
        delta.update(-1.0,r,a);
        delta.scale(1.0/b);

        // y_{n+1} = y_{n} + delta
        y.update(1.0,delta,1.0);

        // Compute residual r = x - A*y
        b_A->apply(y,r);
        r.update(1.0,x,-1.0);

        r.norm2(r_norm());

        if( b_verbosity >= HIGH )
        {
            std::cout << "Relative residual norm at iteration " << b_num_iters
                << " is " << r_norm[0]/r0_norm << std::endl;
        }

        // Check for convergence
        if( r_norm[0]/r0_norm < b_tolerance )
        {
            if( b_verbosity >= LOW )
            {
                std::cout << "Chebyshev Iteration converged after "
                    << b_num_iters << " iterations." << std::endl;
            }
            break;
        }

        // Check for max iterations
        if( b_num_iters >= b_max_iterations )
        {
            if( b_verbosity >= LOW )
            {
                std::cout << "Chebyshev Iteration reached maximum iteration "
                    << "count with relative residual norm of "
                    << r_norm[0]/r0_norm << std::endl;
            }
            break;
        }
    }
}

} // namespace alea

