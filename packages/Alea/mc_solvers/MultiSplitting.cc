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

MultiSplitting::MultiSplitting(Teuchos::RCP<const MATRIX> A,
           Teuchos::RCP<Teuchos::ParameterList> pl )
  : AleaSolver(A,pl)
{
    // Get MultiSplitting pl
    Teuchos::RCP<Teuchos::ParameterList> solver_pl =
        Teuchos::sublist(pl,"MultiSplitting");
        
             
    
    
    // Get InnerSolver pl
    Teuchos::sublist(inner_solver_pl,"InnerSolver");
    
    if(d_inner_solver == "monte_carlo")
    	subdomains.resize( d_num_blocks, MonteCarloSolver );
    	for(unsigned int p = 0; p!=d_num_blocks; ++p)
    	{
    		subdomains[p].part[0]=
    	}
    else
    	subdomains.resize( d_num_blocks, RichardsonIteration );
    	for(unsigned int p = 0; p!=d_num_blocks; ++p)
    	{
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
void MultiSplitting::applyImpl(const MV &x, MV &y) const
{
}

} // namespace alea

