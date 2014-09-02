//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/LinearSolverBuilder.t.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 12:20:24 2014
 * \brief  LinearSolverBuilder template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_LinearSolverBuilder_t_hh
#define solvers_LinearSolverBuilder_t_hh

#include <string>

#include "utils/String_Functions.hh"
#include "LinearSolverBuilder.hh"

#include "StratimikosSolver.hh"
#include "Richardson.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Build a denovo LinearSolver.
 *
 * This function creates a LinearSolver object from a given database.  The
 * logic of selecting a particular solver is designed to maintain backwards
 * compatibility with previous functionality.  First we look for a database
 * entry "solver_type", which can be "profugus" or "stratimikos".  If that
 * entry exists, the corresponding solver type will be built.  If not, we look
 * for database entries "profugus_solver" and build the appropriate class.
 * Current valid "profugus_solver" options are "Richardson".
 *
 */
//---------------------------------------------------------------------------//

template <class MV, class OP>
Teuchos::RCP<LinearSolver<MV,OP> >
LinearSolverBuilder<MV,OP>::build_solver( RCP_ParameterList db )
{
    using std::string;

    RCP_LinearSolver solver;

    // Determine type of solver to be constructed (defaults to profugus)
    string solver_type = to_lower(
        db->get<string>("solver_type", string("profugus")));

    // Check for native solvers
    if (solver_type == "profugus")
    {
        // get profugus solver type
        string type = to_lower(
            db->get<string>("profugus_solver", string("richardson")));

        if (type == "richardson")
        {
            solver = Teuchos::rcp( new Richardson<MV,OP>(db));
        }
        else
        {
            VALIDATE(false, "Invalid 'profugus_solver' type of "
                      << type << " entered.  Valid entries are 'richardson'");
        }
    }
    else if (solver_type == "stratimikos")
    {
        // Just build the stratimikos solver, let validation be handled there
        solver = Teuchos::rcp(new StratimikosSolver<MV,OP>(db));
    }
    else
    {
        VALIDATE(false, "Error: Invalid LinearSolver option "
                 "'" << solver_type << "'\n"
                 "Specify linear solver by setting solver_type="
                 "'profugus' or 'stratimikos' or by setting the "
                 "profugus_solver database entry.\n");
    }

    return solver;
}

} // end namespace profugus

#endif // solvers_LinearSolverBuilder_t_hh

//---------------------------------------------------------------------------//
//                 end of LinearSolverBuilder.t.hh
//---------------------------------------------------------------------------//
