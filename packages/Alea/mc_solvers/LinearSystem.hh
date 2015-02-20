//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSystem.hh
 * \author Steven Hamilton
 * \brief  LinearSystem class declarations.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_LinearSystem_hh
#define mc_solvers_LinearSystem_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

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
        TEUCHOS_ASSERT( d_A != Teuchos::null );
        return d_A;
    }

    //! Return problem right hand side vector
    Teuchos::RCP<const MV> getRhs() const
    {
        TEUCHOS_ASSERT( d_b != Teuchos::null );
        return d_b;
    }

  private:

    // Problem data
    Teuchos::RCP<const MATRIX> d_A;
    Teuchos::RCP<const MV>     d_b;
};

}

#endif // mc_solvers_LinearSystem_hh

