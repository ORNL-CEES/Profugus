//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSystem.hh
 * \author Massimiliano Lupo Pasini
 * \brief  EigenSystem class declarations.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_EigenSystem_hh
#define mc_solvers_EigenSystem_hh

#include "AleaTypedefs.hh"
#include "harness/DBC.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class EigenSystem
 * \brief Container for data to solve eigen problems.
 */
//---------------------------------------------------------------------------//
class EigenSystem
{
  public:

    EigenSystem(Teuchos::RCP<const MATRIX> A,
                 Teuchos::RCP<const MV>     b);

    //! Return problem matrix
    Teuchos::RCP<const MATRIX> getMatrix() const
    {
        REQUIRE( d_A != Teuchos::null );
        return d_A;
    }

    //! Return problem right hand side vector
    Teuchos::RCP<const MV> getInitialGuess() const
    {
        REQUIRE( d_initial_guess != Teuchos::null );
        return d_initial_guess;
    }

  private:

    // Problem data
    Teuchos::RCP<const MATRIX> d_A;
    Teuchos::RCP<const MV>     d_initial_guess;
};

}

#endif // mc_solvers_LinearSystem_hh

