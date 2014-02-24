//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Davidson_Eigensolver.cc
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 14:41:41 2014
 * \brief  Davidson_Eigensolver member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "AnasaziBasicEigenproblem.hpp"

#include "harness/Soft_Equivalence.hh"
#include "comm/P_Stream.hh"
#include "Davidson_Eigensolver.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//

Davidson_Eigensolver::Davidson_Eigensolver(RCP_ParameterList db,
                                           RCP_OP            LHS,
                                           RCP_OP            RHS)
    : Base(db)
    , d_db(db)
    , d_LHS(LHS)
    , d_RHS(RHS)
{
    Require(!d_db.is_null());
    Require(d_LHS!=Teuchos::null);
    Require(d_RHS!=Teuchos::null);

    // make a default Anasazi database
    RCP_ParameterList anasazi_db = Teuchos::sublist(db, "Anasazi");

    // Set desired eigenvalue type
    // We're actually going to solve the eigenproblem as
    //  B*x = k*A*x, because that way we're searching for k
    //  rather than 1/k
    anasazi_db->set("Which", std::string("LM"));

    // Set some default database entries
    anasazi_db->get("Convergence Tolerance",1e-6);
    anasazi_db->get("Maximum Subspace Dimension",25);
    anasazi_db->get("Restart Dimension",5);
    anasazi_db->get("Maximum Restarts",100);
    anasazi_db->get("Initial Guess",std::string("User"));

    // Set verbosity of solver
    anasazi_db->get("Output Level", std::string("low"));
    std::string output_level = anasazi_db->get<std::string>("Output Level");

    // Map Denovo "Output Level" to Anasazi "Verbosity"
    int verbosity = Anasazi::Errors;
    if( output_level == "low" )
        verbosity = Anasazi::FinalSummary;
    else if( output_level == "medium" )
        verbosity = Anasazi::IterationDetails;
    else if( output_level == "high" )
        verbosity = Anasazi::Debug;

    anasazi_db->get("Verbosity", verbosity);
}

//---------------------------------------------------------------------------//
// SOLVE EIGENPROBLEM
//---------------------------------------------------------------------------//

void Davidson_Eigensolver::solve( double &keff, Teuchos::RCP<MV> x)
{
    Require (d_db->isSublist("Anasazi"));

    // Create eigenproblem
    RCP<Anasazi::BasicEigenproblem<double,MV,OP> > problem(
        new Anasazi::BasicEigenproblem<double,MV,OP>() );

    // Remember, we're switching usual convention of LHS and RHS
    //  so that we're converging on k rather than 1/k
    problem->setA(d_RHS);
    problem->setM(d_LHS);
    problem->setPrec(d_prec);
    problem->setInitVec(x);
    problem->setNEV(1);
    bool problem_set = problem->setProblem();
    Ensure( problem_set );

    // Extract Anasazi DB
    ParameterList &anasazi_list = d_db->sublist("Anasazi");

    // Create solver
    Anasazi::GeneralizedDavidsonSolMgr<double,MV,OP> solver(
            problem, anasazi_list);

    // Solve
    Anasazi::ReturnType davReturn = solver.solve();

    Insist( davReturn == Anasazi::Converged,
            "Davidson Solver did not converge!" );

    // Extract solution
    Anasazi::Eigensolution<double,MV> solution =
        solver.getProblem().getSolution();
    Anasazi::Value<double> eval = (solution.Evals)[0];
    Check( profugus::soft_equiv( eval.imagpart, 0.0 ) );
    keff = eval.realpart;

    // Get view of dominant eigenvector
    x->Update(1.0,*(solution.Evecs),0.0);

    if( b_verbosity >= LOW )
    {
        profugus::pout << "+++ Block Generalized Davidson Eigensolver "
                       << "converged in "
                       << solver.getNumIters()
                       << " iterations" << profugus::endl;
    }
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Davidson_Eigensolver.cc
//---------------------------------------------------------------------------//
