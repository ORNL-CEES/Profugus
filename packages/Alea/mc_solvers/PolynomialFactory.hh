
#ifndef Alea_mc_solvers_PolynomialFactory_hh
#define Alea_mc_solvers_PolynomialFactory_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "AnasaziTypes.hpp"

#include "AleaTypedefs.hh"

#include "Polynomial.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class PolynomialFactory
 * \brief Create Polynomial object
 */
//---------------------------------------------------------------------------//
class PolynomialFactory
{

  private:

    // No construction
    PolynomialFactory(){};

  public:

    static Teuchos::RCP<Polynomial>
    buildPolynomial(Teuchos::RCP<const MATRIX> A,
                    Teuchos::RCP<Teuchos::ParameterList> pl );
};

}

#endif // Alea_mc_solvers_PolynomialFactory_hh

