//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSystem.cc
 * \author Steven Hamilton
 * \brief  LinearSystem class definitions.
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "LinearSystem.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param A  Problem matrix
 * \param b  Problem right hand side vector
 */
//---------------------------------------------------------------------------//
LinearSystem::LinearSystem(Teuchos::RCP<const MATRIX> A,
                           Teuchos::RCP<const MV>     b)
  : d_A(A)
  , d_b(b)
{
}

} // namespace alea

