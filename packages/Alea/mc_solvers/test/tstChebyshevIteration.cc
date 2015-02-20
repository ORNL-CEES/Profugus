//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   testChebyshevIteration.cc
 * \author Steven Hamilton
 * \brief  Test of ChebyshevIteration class.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <time.h>

#include "comm/global.hh"
#include "../LinearSystemFactory.hh"
#include "../ChebyshevIteration.hh"
#include "../AleaTypedefs.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"

using namespace alea;

TEST(Chebyshev, Basic)
{
    // Read ParameterList from file
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(pl,"Problem");

    mat_pl->set("matrix_type","laplacian");
    mat_pl->set("matrix_size",20);
    pl->set("max_iterations",1000);
    pl->set("tolerance",1.0e-6);
    pl->sublist("Chebyshev").set("verbosity","high");

    Teuchos::RCP<LinearSystem> system =
        alea::LinearSystemFactory::buildLinearSystem(pl);
    Teuchos::RCP<const MATRIX> A = system->getMatrix();
    Teuchos::RCP<const MV> b = system->getRhs();

    // Run Chebyshev solver
    Teuchos::RCP<alea::ChebyshevIteration> solver(
        new alea::ChebyshevIteration(A,pl) );

    Teuchos::RCP<MV> x( new MV(A->getDomainMap(),1) );
    x->putScalar(0.0);

    solver->apply(*b,*x);

    LO num_iters = solver->getNumIters();
    EXPECT_EQ( 110, num_iters );

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

