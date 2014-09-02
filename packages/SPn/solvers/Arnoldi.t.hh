//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Arnoldi.cc
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 14:41:30 2014
 * \brief  Arnoldi template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_Arnoldi_t_hh
#define solvers_Arnoldi_t_hh

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
template <class MV, class OP> Arnoldi<MV,OP>::Arnoldi(RCP_ParameterList db)
    : EigenvalueSolver<MV,OP>(db)
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
    ENSURE(!d_pl.is_null());
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 *\brief Set operator for eigensolver.
 */
template <class MV, class OP>
void Arnoldi<MV,OP>::set_operator(RCP_OP A)
{
    REQUIRE(!A.is_null());
    d_A = A;

}

//---------------------------------------------------------------------------//
/*!
 *\brief Solve eigenproblem.
 *
 * This function solves the eigenproblem defined by set_operator.  The
 * dominant eigenvalue will be placed in the input argument eval and the
 * corresponding eigenvector in evec.
 */
template <class MV, class OP>
void Arnoldi<MV,OP>::solve(double &eval,
                    RCP_MV  evec)
{
    REQUIRE(!evec.is_null());
    REQUIRE(!d_A.is_null());

    // Create eigenproblem
    RCP_Eigenproblem problem(new Eigenproblem(d_A, evec));
    problem->setNEV(1);
    problem->setProblem();

    // Set the appropriate number of blocks based on requested subspace
    if( !d_pl->isType<int>("Num Blocks") )
    {
        // calculate problem sizes
        int subspace   = d_pl->get<int>("subspace");
        int N = MultiVecTraits::GetVecLength(*evec);
        int block_size = d_pl->get<int>("Block Size");
        int num_blocks = std::min( (N-2)/block_size, subspace );

        // add remaining defaults
        d_pl->set("Num Blocks", num_blocks);
    }

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
    CHECK( MultiVecTraits::GetNumberVecs(*outvec) > 0 );

    // Assign the first vector of the eigensolution (outvec may contain
    //  several even though we only converged on one) to the first vector
    //  of evec.
    std::vector<int> ind(1,0);
    MultiVecTraits::SetBlock(*outvec,ind,*evec);
}

} // end namespace profugus

#endif // solvers_Arnoldi_t_hh

//---------------------------------------------------------------------------//
//                 end of Arnoldi.t.hh
//---------------------------------------------------------------------------//
