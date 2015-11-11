//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MultiSplitting.cc
 * \author Massimiliano Lupo Pasini
 * \brief  Perform Multi Splitting with either deterministic Richardson 
 * iterations or MC updates.
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "Teuchos_ArrayRCP.hpp"

#include "MultiSplitting.hh"
#include "LinearSystem_MultiSplitting.hh"
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
 * "MultiSplitting"
 * sublist:
 *  - block_size(int)             : >0 (1)
 *  - overlapping(SCALAR)         : >0.0 (0.0)
 */
//---------------------------------------------------------------------------//

MultiSplitting::MultiSplitting( Teuchos::RCP<Teuchos::ParameterList> &pl )

{
    // Get MultiSplitting pl
    d_pl = pl;
    Teuchos::RCP<LinearSystem_MultiSplitting> ms( new LinearSystem_MultiSplitting(d_pl) );
    d_multisplitting = ms;                
             
    d_A = d_multisplitting->getMatrix();
    std::cout<< d_A->getGlobalNumRows()<<std::endl;
    std::cout<<"Is upper triangular: "<<d_A->isUpperTriangular()<<std::endl;
    std::cout<<"Is lower triangular: "<<d_A->isLowerTriangular()<<std::endl;

    d_b = d_multisplitting->getRhs();         
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(d_pl,"MultiSplitting");             
             
    d_inner_solver = d_multisplitting->getInnerSolverType();
    VALIDATE( d_inner_solver == "richardson" || d_inner_solver == "monte_carlo",
         "the only iterative solvers admitted are Richardson and MonteCarlo" );         
             
    d_divergence_tol = mat_pl->get("divergence_tolerance",1.0e4);
    std::cout<<"divergence tolerance "<<d_divergence_tol<<std::endl;
    
    b_max_iterations = mat_pl->get<int>("max_iterations");
    
    if( mat_pl->isType<std::string>("verbosity") )
    {
        std::string verbosity = pl->get<std::string>("verbosity");
        VALIDATE(verbosity=="none"   || verbosity=="low"  ||
                 verbosity=="medium" || verbosity=="high" ||
                 verbosity=="debug",
                 "Invalid verbosity specified.");
        if( verbosity == "none")
        {
            b_verbosity = NONE;
        }
        else if( verbosity == "low" )
        {
            b_verbosity = LOW;
        }
        else if( verbosity == "medium" )
        {
            b_verbosity = MEDIUM;
        }
        else if( verbosity == "high" )
        {
            b_verbosity = HIGH;
        }
        else if( verbosity == "debug" )
        {
            b_verbosity = DEBUG;
        }
    }
    
    
}


//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Perform MultiSplitting.
 */
//---------------------------------------------------------------------------//
void MultiSplitting::solve(Teuchos::RCP<MV> &x) const
{
    splitting split = d_multisplitting->buildSplitting(d_pl,0);
 
    Teuchos::RCP<alea::AleaSolver> solver =
        alea::LinearSolverFactory::buildSolver(d_inner_solver,
         split.A,d_pl);

    solver->apply(*(split.b),*x);

    // For now we only support operating on a single vector
    REQUIRE( x->getNumVectors() == 1 );
    // Compute initial residual
    MV r( (split.b)->getMap(),1 );

    Teuchos::ArrayRCP<MAGNITUDE> r_norm(1);
    r.norm2(r_norm());
    if( r_norm[0] == 0.0 )
        return;
        
        MAGNITUDE r0_norm = r_norm[0];

    b_num_iters = 0;
    while( true )
    {
        b_num_iters++;

        // Compute residual r = x - A*y
        d_A->apply(*x,r);
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
            std::cout << "Richardson Iteration reached maximum iteration "
                 << "count with relative residual norm of "
                 << r_norm[0]/r0_norm << std::endl;
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

