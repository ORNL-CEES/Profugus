//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/LinearSystem.hh
 * \author Steven Hamilton
 * \brief  LinearSystem class declarations.
 */
//---------------------------------------------------------------------------//

#ifndef Alea_mc_solvers_LinearSystem_hh
#define Alea_mc_solvers_LinearSystem_hh

#include "AleaTypedefs.hh"
#include "harness/DBC.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class LinearSystem
 * \brief Container for data to solve linear system.
 */
//---------------------------------------------------------------------------//
class LinearSystem
{
  public:

    LinearSystem(Teuchos::RCP<const MATRIX> A,
                 Teuchos::RCP<const MV>     b);

    //! Return problem matrix
    Teuchos::RCP<const MATRIX> getMatrix() const
    {
        REQUIRE( d_A != Teuchos::null );
        return d_A;
    }

    //! Return problem right hand side vector
    Teuchos::RCP<const MV> getRhs() const
    {
        REQUIRE( d_b != Teuchos::null );
        return d_b;
    }

  private:

    // Problem data
    Teuchos::RCP<const MATRIX> d_A;
    Teuchos::RCP<const MV>     d_b;
};

}

#endif // Alea_mc_solvers_LinearSystem_hh

