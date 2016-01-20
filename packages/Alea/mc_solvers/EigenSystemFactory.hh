//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSystemFactory.hh
 * \author Massimiliano Lupo Pasini
 * \brief  Construct EpetraCrsMatrix from ParameterList.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_EigenSystemFactory_hh
#define mc_solvers_EigenSystemFactory_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "EigenSystem.hh"
#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class EigenSystemFactory
 * \brief Load/construct matrix.
 *
 * This class provides a simple interface to build an eigenvalue problem. Several
 * simple matrices can be constructed from the provided ParameterList or
 * a matrix can be loaded from a Matrix Market file.  An initial guess for
 * the linear system will also be constructed, by assigning a dummy vector if no
 * meaningful vector is available or assuming n initial guess from a file.
 */
//---------------------------------------------------------------------------//
class EigenSystemFactory
{
  private:

    // Pure static, disallow construction
    EigenSystemFactory(){};

  public:

    static Teuchos::RCP<EigenSystem> buildEigenSystem(
        Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    static void buildDiffusionSystem(
        Teuchos::RCP<Teuchos::ParameterList> pl,
        Teuchos::RCP<CRS_MATRIX> &A,
        Teuchos::RCP<MV>     &b );

    static void buildConvectionDiffusionSystem(
        Teuchos::RCP<Teuchos::ParameterList> pl,
        Teuchos::RCP<CRS_MATRIX> &A,
        Teuchos::RCP<MV>     &b );

    static Teuchos::RCP<CRS_MATRIX> buildLaplacianMatrix(
        int N, Teuchos::RCP<Teuchos::ParameterList> pl);

    static Teuchos::RCP<CRS_MATRIX> buildConvectionMatrix(
        int N, Teuchos::RCP<Teuchos::ParameterList> pl);

    static void buildMatrixMarketSystem(
        Teuchos::RCP<Teuchos::ParameterList> pl,
        Teuchos::RCP<CRS_MATRIX> &A,
        Teuchos::RCP<MV>     &b );

    static void buildProfugusSystem(
        Teuchos::RCP<Teuchos::ParameterList> pl,
        Teuchos::RCP<CRS_MATRIX> &A,
        Teuchos::RCP<MV>     &b );

    static Teuchos::RCP<CRS_MATRIX> applyShift(
        Teuchos::RCP<CRS_MATRIX>                 A,
        Teuchos::RCP<Teuchos::ParameterList> pl);

};

}

#endif // mc_solvers_EigenSystemFactory_hh

