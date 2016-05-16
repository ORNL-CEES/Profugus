//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/SpnSolverBuilder.cc
 * \author Thomas M. Evans
 * \date   Tue Oct 23 21:17:05 2012
 * \brief  SpnSolverBuilder member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "harness/DBC.hh"
#include "SpnSolverBuilder.hh"

#include "solvers/LinAlgTypedefs.hh"
#include "Eigenvalue_Solver.hh"
#include "Fixed_Source_Solver.hh"
#include "Time_Dependent_Solver.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * The number of equations is equal to \f$(N+1)/2\f$.
 *
 * \param N SPN order (1, 3, 5, 7)
 */
Teuchos::RCP<Solver_Base> SpnSolverBuilder::build(std::string       problem,
                                                  RCP_ParameterList db)
{
    REQUIRE( (problem=="eigenvalue") || (problem=="fixed") ||
             (problem=="fixed_tdep") );

    // Determine if epetra or tpetra should be used
    std::string implementation = db->get("trilinos_implementation",
        std::string("epetra"));
    REQUIRE( implementation == "epetra" || implementation == "tpetra" );

    Teuchos::RCP<Solver_Base> solver;

    if( problem == "eigenvalue" )
    {
        if( implementation == "epetra" )
            solver = Teuchos::rcp( new Eigenvalue_Solver<EpetraTypes>(db) );
        else if( implementation == "tpetra" )
            solver = Teuchos::rcp( new Eigenvalue_Solver<TpetraTypes>(db) );
    }
    else if( problem == "fixed" )
    {
        if( implementation == "epetra" )
            solver = Teuchos::rcp( new Fixed_Source_Solver<EpetraTypes>(db) );
        else if( implementation == "tpetra" )
            solver = Teuchos::rcp( new Fixed_Source_Solver<TpetraTypes>(db) );
    }
    else if( problem == "fixed_tdep" )
    {
        if( implementation == "epetra" )
            solver = Teuchos::rcp( new Time_Dependent_Solver<EpetraTypes>(db) );
        else if( implementation == "tpetra" )
            solver = Teuchos::rcp( new Time_Dependent_Solver<TpetraTypes>(db) );
    }

    ENSURE( solver != Teuchos::null );
    return solver;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of SpnSolverBuilder.cc
//---------------------------------------------------------------------------//
