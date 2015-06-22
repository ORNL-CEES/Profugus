//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SeuqntialMC.cc
 * \author Steven Hamilton
 * \brief  Performs Sequential Monte Carlo on linear system
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "Teuchos_ArrayRCP.hpp"

#include "SequentialMC.hh"
#include "LinearSolverFactory.hh"
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
 * Behavior is controlled by following PL entries on the nested
 * "SequentialMC"
 * sublist:
 *  - max_iterations(int)             : >0 (1000)
 *  - damping(SCALAR)                 : >0.0 (1.0)
 *  - tolerance(MAGNITUDE)            : >0.0 (1.0e-6)
 *  - divergence_tolerance(MAGNITUDE) : >0.0 (1.0e4),
 *                                      residual norm to declare failure
 *  - verbosity(string)               : "none", ("low"), "medium", "high"
 */
//---------------------------------------------------------------------------//
SequentialMC::SequentialMC(Teuchos::RCP<const MATRIX> A,
           Teuchos::RCP<Teuchos::ParameterList> pl )
  : AleaSolver(A,pl)
{
    // Get SequentialMC pl
    Teuchos::RCP<Teuchos::ParameterList> solver_pl =
        Teuchos::sublist(pl,"SequentialMC");

    d_damping = solver_pl->get("damping_factor",1.0);

    d_divergence_tol = solver_pl->get("divergence_tolerance",1.0e4);

    // Override settings if present on local list
    this->setParameters(solver_pl);

    // Build preconditioner, default to monte_carlo (MCSA)
    std::string prec_type = "monte_carlo"; // there is no point in admitting other preconditioners but MC in this case
    d_preconditioner = LinearSolverFactory::buildSolver(prec_type,A,pl);

    b_label = "SequentialMC";
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Perform SyntheticAcceleration iterations.
 */
//---------------------------------------------------------------------------//
void SequentialMC::applyImpl(const MV &x, MV &y) const
{
    // For now we only support operating on a single vector
    REQUIRE( x.getNumVectors() == 1 );

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
    
    y.update(1.0,r,1.0);
    b_A->apply(y,r);
    r.update(1.0,x,-1.0);

    while( true )
    {
      	b_num_iters++;

        Pr.update(1.0,r,0.0);
        d_preconditioner->apply(r,Pr);

        // Update y: y = y + Pr
        y.update(d_damping,Pr,1.0);

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
                std::cout << "SequentialMC converged after "
                    << b_num_iters << " iterations." << r_norm[0]/r0_norm << std::endl;
            }
            break;
        }

        // Check for max iterations
        if( b_num_iters >= b_max_iterations )
        {
            if( b_verbosity >= LOW )
            {
                std::cout << "SyntheticAcceleration reached maximum iteration "
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
                std::cout << "SyntheticAcceleration diverged after "
                    << b_num_iters << " iterations." << std::endl;
            }
            b_num_iters = -b_num_iters;
            break;
        }
    }
}

} // namespace alea

