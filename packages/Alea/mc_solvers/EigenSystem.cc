//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSystem.cc
 * \author Massimiliano Lupo Pasini
 * \brief  EigenSystem class definitions.
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "EigenSystem.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param A  Problem matrix
 * \param b  Problem initial guess vector
 */
//---------------------------------------------------------------------------//
EigenSystem::EigenSystem(Teuchos::RCP<const MATRIX> A,
                           Teuchos::RCP<const MV>     b)
  : d_A(A)
  , d_initial_guess(b)
{
}

} // namespace alea

