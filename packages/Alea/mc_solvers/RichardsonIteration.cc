//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RichardsonIteration.cc
 * \author Steven Hamilton
 * \brief  Perform Adjoint Monte Carlo on linear system
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "RichardsonIteration.hh"
#include "LinearSolverFactory.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param A  Problem matrix
 * \param pl ParameterList
 *
 * Behavior is controlled by following PL entries on the nested "Richardson"
 * sublist:
 *  - preconditioner(string) : ("none") any valid AleaSolver
 *  - max_iterations(int)    : >0 (1000)
 *  - tolerance(SCALAR)      : >0.0 (1.0e-6)
 *  - verbosity(string)      : "none", ("low"), "medium", "high"
 */
//---------------------------------------------------------------------------//
RichardsonIteration::RichardsonIteration(
        Teuchos::RCP<const MATRIX> A,
        Teuchos::RCP<Teuchos::ParameterList> pl )
  : AleaSolver(A,pl)
{
    // Get belos pl
    Teuchos::RCP<Teuchos::ParameterList> rich_pl =
        Teuchos::sublist(pl,"Richardson");

    d_divergence_tol = rich_pl->get("divergence_tolerance",1.0e4);

    // Override default parameters if present on sublist
    this->setParameters(rich_pl);

    std::string prec_type = pl->get("preconditioner","none");

    d_P = LinearSolverFactory::buildSolver(prec_type,A,pl);

    b_label = "RichardsonIteration";
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Perform Richardson iteration.
 */
//---------------------------------------------------------------------------//
void RichardsonIteration::applyImpl(const MV &x, MV &y) const
{
    // For now we only support operating on a single vector
    TEUCHOS_TEST_FOR_EXCEPT( x.getNumVectors() != 1 );

    // Compute initial residual
    MV r(y.getMap(),1);
    r.update(1.0,x,0.0);
    MV Pr(y.getMap(),1);

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

    b_num_iters = 0;
    while( true )
    {
        b_num_iters++;

        // Precondition residual
        if( d_P != Teuchos::null )
        {
            d_P->apply(r,Pr,Teuchos::NO_TRANS,1.0,0.0);
        }
        else
        {
            Pr.update(1.0,r,0.0);
        }

        // y = y + Pr
        y.update(1.0,Pr,1.0);

        // Compute residual r = x - A*y
        b_A->apply(y,r);
        r.update(1.0,x,-1.0);

        // Check convergence on true (rather than preconditioned) residual
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
                std::cout << "Richardson Iteration converged after "
                    << b_num_iters << " iterations." << std::endl;
            }
            break;
        }

        // Check for max iterations
        if( b_num_iters >= b_max_iterations )
        {
            if( b_verbosity >= LOW )
            {
                std::cout << "Richardson Iteration reached maximum iteration "
                    << "count with relative residual norm of "
                    << r_norm[0]/r0_norm << std::endl;
            }
            b_num_iters = -1;
            break;
        }

        // Check for divergence
        if( r_norm[0]/r0_norm > d_divergence_tol)
        {
            if( b_verbosity >= LOW )
            {
                std::cout << "Richardson Iteration diverged after "
                    << b_num_iters << " iterations." << std::endl;
            }
            b_num_iters = -b_num_iters;
            break;
        }

    }
}

} // namespace alea

