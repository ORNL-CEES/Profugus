//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   testRichardsonIteration.cc
 * \author Steven Hamilton
 * \brief  Test of RichardsonIteration class.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <time.h>

#include "comm/global.hh"
#include "../LinearSystemFactory.hh"
#include "../RichardsonIteration.hh"
#include "../AleaTypedefs.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"

using namespace alea;

TEST(Richardson, Basic)
{
    // Read ParameterList from file
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(pl,"Problem");

    mat_pl->set("matrix_type","laplacian");
    mat_pl->set("matrix_size",10);
    pl->set("max_iterations",1000);
    pl->set("tolerance",1.0e-6);
    pl->sublist("Richardson").set("verbosity","high");

    Teuchos::RCP<LinearSystem> system =
        alea::LinearSystemFactory::buildLinearSystem(pl);
    Teuchos::RCP<const MATRIX> A = system->getMatrix();
    Teuchos::RCP<const MV> b = system->getRhs();

    // Test collision estimator
    Teuchos::RCP<alea::RichardsonIteration> solver(
        new alea::RichardsonIteration(A,pl) );

    Teuchos::RCP<MV> x( new MV(A->getDomainMap(),1) );
    x->putScalar(0.0);

    solver->apply(*b,*x);

    LO num_iters = solver->getNumIters();
    EXPECT_EQ( 277, num_iters );

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

    // Reset and use Neumann polynomial preconditioning
    pl->set("preconditioner","polynomial");
    Teuchos::RCP<Teuchos::ParameterList> poly_pl =
        Teuchos::sublist(pl,"Polynomial");
    poly_pl->set("polynomial_type","neumann");
    poly_pl->set("polynomial_order",4);
    solver = Teuchos::rcp( new alea::RichardsonIteration(A,pl) );

    // Solve
    x->putScalar(0.0);
    solver->apply(*b,*x);

    num_iters = solver->getNumIters();

    EXPECT_EQ( 56, num_iters );

    A->apply(*x,*r);
    r->update(1.0,*b,-1.0);
    r->norm2(res_norm());
    std::cout << "Final relative residual norm: "
              << res_norm[0]/b_norm[0] << std::endl;
    EXPECT_TRUE( res_norm[0]/b_norm[0] < 1.0e-6 );

    // One more time with Chebyshev polynomial preconditioning
    poly_pl->set("polynomial_type","chebyshev");
    solver = Teuchos::rcp( new alea::RichardsonIteration(A,pl) );

    // Solve
    x->putScalar(0.0);
    solver->apply(*b,*x);

    num_iters = solver->getNumIters();

    EXPECT_EQ( 15, num_iters );

    A->apply(*x,*r);
    r->update(1.0,*b,-1.0);
    r->norm2(res_norm());
    std::cout << "Final relative residual norm: "
              << res_norm[0]/b_norm[0] << std::endl;
    EXPECT_TRUE( res_norm[0]/b_norm[0] < 1.0e-6 );
}

