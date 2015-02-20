//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   testMonteCarloSolver.cc
 * \author Steven Hamilton
 * \brief  Test of MonteCarloSolver class.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <time.h>

#include "../LinearSystem.hh"
#include "../LinearSystemFactory.hh"
#include "../MonteCarloSolver.hh"
#include "../DeviceTraits.hh"
#include "../AleaTypedefs.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"

using namespace alea;

class MonteCarlo : public ::testing::Test
{
  protected:

    void SetUp()
    {
        DeviceTraits<DEVICE>::initializeDevice(1);

        // Create ParameterList
        d_pl = Teuchos::rcp( new Teuchos::ParameterList() );
        Teuchos::RCP<Teuchos::ParameterList> mat_pl =
            Teuchos::sublist(d_pl,"Problem");
        Teuchos::RCP<Teuchos::ParameterList> mc_pl =
            Teuchos::sublist(d_pl,"Monte Carlo");
        Teuchos::RCP<Teuchos::ParameterList> poly_pl =
            Teuchos::sublist(d_pl,"Polynomial");

        mat_pl->set("matrix_type","laplacian");
        mat_pl->set("matrix_size",10);

        mc_pl->set("num_histories",20000);
        mc_pl->set("weight_cutoff",0.01);
        mc_pl->set("verbosity","medium");

        poly_pl->set("polynomial_order",100);

        Teuchos::RCP<alea::LinearSystem> system =
            alea::LinearSystemFactory::buildLinearSystem(d_pl);
        d_A = system->getMatrix();
        d_b = system->getRhs();
    }

    void TearDown()
    {
        DeviceTraits<DEVICE>::finalizeDevice();
    }

    void Solve(double expected_tol)
    {
        Teuchos::RCP<alea::MonteCarloSolver> solver(
            new alea::MonteCarloSolver(d_A,d_pl) );
        solver->compute();

        Teuchos::RCP<MV> x( new MV(d_A->getDomainMap(),1) );

        solver->apply(*d_b,*x);

        // Compute final residual
        Teuchos::RCP<MV> r( new MV(d_A->getDomainMap(),1) );
        d_A->apply(*x,*r);
        r->update(1.0,*d_b,-1.0);
        Teuchos::ArrayRCP<SCALAR> res_norm(1), b_norm(1);
        r->norm2(res_norm());
        d_b->norm2(b_norm());
        std::cout << "Final relative residual norm: "
                  << res_norm[0]/b_norm[0] << std::endl;

        // This should *almost* always pass
        EXPECT_TRUE( res_norm[0]/b_norm[0] < expected_tol );
    }

    Teuchos::RCP<Teuchos::ParameterList> d_pl;
    Teuchos::RCP<const MATRIX> d_A;
    Teuchos::RCP<const MV> d_b;
    Teuchos::RCP<MV> d_x;
};

TEST_F(MonteCarlo,Collision)
{
    Teuchos::sublist(d_pl,"Monte Carlo")->set("estimator","collision");
    this->Solve(0.12);
}

TEST_F(MonteCarlo,ExpectedValue)
{
    Teuchos::sublist(d_pl,"Monte Carlo")->set("estimator","expected_value");
    this->Solve(0.06);
}
