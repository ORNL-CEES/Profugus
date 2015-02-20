//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   testMCSA.cc
 * \author Steven Hamilton
 * \brief  Test of MCSA class.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <time.h>

#include "../LinearSystemFactory.hh"
#include "../SyntheticAcceleration.hh"
#include "../AleaTypedefs.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"

using namespace alea;

TEST(MCSA, Basic)
{
    // Read ParameterList from file
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(pl,"Problem");
    Teuchos::RCP<Teuchos::ParameterList> mc_pl =
        Teuchos::sublist(pl,"Monte Carlo");
    Teuchos::RCP<Teuchos::ParameterList> poly_pl =
        Teuchos::sublist(pl,"Polynomial");

    mat_pl->set("matrix_type","laplacian");
    mat_pl->set("matrix_size",25);

    mc_pl->set("num_histories",1000);
    mc_pl->set("weight_cutoff",0.1);
    mc_pl->set("estimator","expected_value");

    poly_pl->set("polynomial_order",20);

    pl->set("max_iterations",1000);
    pl->set("tolerance",1.0e-6);

    Teuchos::RCP<LinearSystem> system =
        alea::LinearSystemFactory::buildLinearSystem(pl);
    Teuchos::RCP<const MATRIX> A = system->getMatrix();
    Teuchos::RCP<const MV> b = system->getRhs();

    // Build MCSA solver
    Teuchos::RCP<alea::SyntheticAcceleration> solver(
        new alea::SyntheticAcceleration(A,pl) );

    Teuchos::RCP<MV> x( new MV(A->getDomainMap(),1) );

    solver->apply(*b,*x);

    // Compute final residual
    Teuchos::RCP<MV> r( new MV(A->getDomainMap(),1) );
    A->apply(*x,*r);
    r->update(1.0,*b,-1.0);
    Teuchos::ArrayRCP<SCALAR> res_norm(1), b_norm(1);
    r->norm2(res_norm());
    b->norm2(b_norm());
    std::cout << "Final relative residual norm: "
              << res_norm[0]/b_norm[0] << std::endl;

    // Should meet convergence criteria
    EXPECT_TRUE( res_norm[0]/b_norm[0] < 1.0e-6 );
}

