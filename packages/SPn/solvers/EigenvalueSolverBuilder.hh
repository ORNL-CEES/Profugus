//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/EigenvalueSolverBuilder.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Mon Feb 24 13:49:22 2014
 * \brief  EigenvalueSolverBuilder class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_EigenvalueSolverBuilder_hh
#define solvers_EigenvalueSolverBuilder_hh

#include "EigenvalueSolver.hh"

#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"

namespace profugus
{

//===========================================================================//
/*!
 * \class EigenvalueSolverBuilder
 * \brief Factory class for creating EigenvalueSolver.
 *
 * This class provides a generic interface for creating an eigensolver
 * involving an Epetra operator (or operators)
 *
 * \sa EigenvalueSolverBuilder.cc for detailed descriptions.
 */
/*!
 * \example solvers/test/tstEigenvalueSolverBuilder.cc
 *
 * Test of EigenvalueSolverBuilder.
 */
//===========================================================================//

class EigenvalueSolverBuilder
{
  public:
    //@{
    //! Typedefs.
    typedef Epetra_MultiVector                    MV;
    typedef Epetra_Operator                       OP;
    typedef EigenvalueSolver<MV,OP>               EigenvalueSolver_t;
    typedef Teuchos::RCP<EigenvalueSolver_t>      RCP_EigenvalueSolver;
    typedef EigenvalueSolver_t::RCP_ParameterList RCP_ParameterList;
    //@}

    // Build standard eigenvalue solver
    static RCP_EigenvalueSolver build_solver(RCP_ParameterList db,
                                             Teuchos::RCP<OP> A);

    // Build generalized eigenvalue solver with preconditioner
    static RCP_EigenvalueSolver build_solver(
            RCP_ParameterList db,
            Teuchos::RCP<OP> A,
            Teuchos::RCP<OP> B,
            Teuchos::RCP<OP> P = Teuchos::null);
};

} // end namespace profugus

#endif // solvers_EigenvalueSolverBuilder_hh

//---------------------------------------------------------------------------//
//                 end of EigenvalueSolverBuilder.hh
//---------------------------------------------------------------------------//
