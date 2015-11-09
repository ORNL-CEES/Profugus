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

MultiSplitting::MultiSplitting( Teuchos::RCP<Teuchos::ParameterList> pl )

{
    // Get MultiSplitting pl
    Teuchos::RCP<Teuchos::ParameterList> b_pl = pl;
    Teuchos::RCP<LinearSystem_MultiSplitting> ms( new LinearSystem_MultiSplitting(pl) );
    d_multisplitting = ms;                
             
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(pl,"MultiSplitting");             
             
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
void MultiSplitting::applyImpl(const MV &x, MV &y) const
{
}

} // namespace alea

