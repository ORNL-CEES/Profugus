//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/GmresPolynomial.hh
 * \author Steven Hamilton
 * \brief  GmresPolynomial class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef Alea_mc_solvers_GmresPolynomial_hh
#define Alea_mc_solvers_GmresPolynomial_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

#include "PolynomialBasis.hh"
#include "Polynomial.hh"

#include "AnasaziTypes.hpp"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class GmresPolynomial
 * \brief Gmres series polynomial.
 */
//---------------------------------------------------------------------------//
class GmresPolynomial : public Polynomial
{
  public:

    // Constructor
    GmresPolynomial(Teuchos::RCP<const MATRIX> A,
               Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    void buildGmresPolyFromArnoldi();
    void buildFomPolyFromArnoldi();
    void buildGmresPolyFromNormalEqns();
    void buildGmresPolyFromQrDecomp();

    Teuchos::RCP<const MV> computeInitialKrylovVector() const;

    Teuchos::ArrayRCP<const Anasazi::Value<SCALAR> >
    computeGmresRoots(
        Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SCALAR> > H) const;

    Teuchos::ArrayRCP<const Anasazi::Value<SCALAR> >
    computeFomRoots(
        Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SCALAR> > H) const;

    Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,SCALAR> >
    computeArnoldiMatrix(int order) const;

    Teuchos::ArrayRCP<SCALAR> computeCoeffsFromRoots(
        Teuchos::ArrayRCP<const Anasazi::Value<SCALAR> > roots) const;

    Teuchos::RCP<MV> computeKrylovStandardBasis(int order,
        Teuchos::RCP<const MV> b) const;

};

}

#endif // Alea_mc_solvers_GmresPolynomial_hh

