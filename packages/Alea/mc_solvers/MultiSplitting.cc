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

}

} // namespace alea

