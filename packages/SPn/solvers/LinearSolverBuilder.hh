//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/LinearSolverBuilder.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 12:20:24 2014
 * \brief  LinearSolverBuilder class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_LinearSolverBuilder_hh
#define solvers_LinearSolverBuilder_hh

#include "LinearSolver.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class LinearSolverBuilder
 * \brief Factory class for creating LinearSolver.
 *
 * This class provides a generic interface for creating a linear system
 * involving an Epetra/Tpetra operator.
 *
 * \sa LinearSolverBuilder.cc for detailed descriptions.
 */
/*!
 * \example solvers/test/tstLinearSolverBuilder.cc
 *
 * Test of LinearSolverBuilder.
 */
//===========================================================================//

template <class MV, class OP>
class LinearSolverBuilder
{
  public:
    //@{
    //! Typedefs.
    typedef LinearSolver<MV,OP>                  LinearSolver_t;
    typedef Teuchos::RCP<LinearSolver_t>         RCP_LinearSolver;
    typedef Teuchos::RCP<Teuchos::ParameterList> RCP_ParameterList;
    //@}

    static RCP_LinearSolver build_solver( RCP_ParameterList db );
};

} // end namespace profugus

#endif // solvers_LinearSolverBuilder_hh

//---------------------------------------------------------------------------//
//                 end of LinearSolverBuilder.hh
//---------------------------------------------------------------------------//
