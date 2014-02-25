//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Arnoldi.cc
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 14:41:30 2014
 * \brief  Arnoldi member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <string>
#include <algorithm>

#include "harness/Warnings.hh"
#include "comm/global.hh"
#include "Arnoldi.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 *\brief Constructor with user defined tolerance and max iterations.
 */
Arnoldi::Arnoldi(RCP_ParameterList db)
    : EigenvalueSolver<Epetra_MultiVector,Epetra_Operator>(db)
{
    b_label = "Arnoldi";

    // make default Anasazi database
    d_pl = Teuchos::sublist(db, "Anasazi");

    // add a subspace default
    d_pl->get("subspace", 20);

    // add default block size
    d_pl->get("Block Size", 1);

    d_pl->get("Which", std::string("LM"));
    d_pl->get("Maximum Restarts", 20);
    d_pl->get("Convergence Tolerance", 1.0e-6);
    d_pl->get("Relative Convergence Tolerance", true);
    d_pl->get("Step Size", 1);
    d_pl->get("Verbosity", static_cast<int>(Anasazi::IterationDetails));

    // Create Teuchos parameter list to give to eigensolver
    Ensure (!d_pl.is_null());
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 *\brief Set operator for eigensolver.
 */
void Arnoldi::set_operator(RCP_OP A)
{
    Require(!A.is_null());
    d_A = A;

    // calculate problem sizes
    int subspace   = d_pl->get<int>("subspace");
    int N          = A->OperatorDomainMap().NumGlobalElements();
    int block_size = d_pl->get<int>("Block Size");
    int num_blocks = std::min( (N-2)/block_size, subspace );

    // add remaining defaults
    if( !d_pl->isType<int>("Num Blocks") )
        d_pl->set("Num Blocks", num_blocks);
}

//---------------------------------------------------------------------------//
/*!
 *\brief Solve eigenproblem.
 *
 * This function solves the eigenproblem defined by set_operator.  The
 * dominant eigenvalue will be placed in the input argument eval and the
 * corresponding eigenvector in evec.
 */
void Arnoldi::solve(double &eval,
                    RCP_MV  evec)
{
    Require (!evec.is_null());
    Require (!d_A.is_null());

    // Create eigenproblem
    RCP_Eigenproblem problem(new Eigenproblem(d_A, evec));
    problem->setNEV(1);
    problem->setProblem();

    // Create eigensolver
    KrylovSchur solver( problem, *d_pl );

    // Solve eigenproblem.
    int returnval = solver.solve();

    // Ensure convergence
    if( returnval == Anasazi::Converged )
    {
        b_converged = true;
    }
    else
    {
        if( profugus::node()==0 )
            ADD_WARNING("Arnoldi failed to converge");
    }

    // Get solution from eigenproblem
    eval          = problem->getSolution().Evals[0].realpart;
    RCP_MV outvec = problem->getSolution().Evecs;
    Check( outvec->NumVectors() > 0 );

    // Assign the first vector of the eigensolution (outvec may contain
    //  several even though we only converged on one) to the first vector
    //  of evec.
    *(*evec)(0) = *(*outvec)(0);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Arnoldi.cc
//---------------------------------------------------------------------------//
