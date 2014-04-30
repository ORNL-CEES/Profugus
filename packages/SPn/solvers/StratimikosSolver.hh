//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/StratimikosSolver.hh
 * \author Chris Baker
 * \date   Mon Sep 17 15:21:00 2012
 * \brief  StratimikosSolver class definition
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_StratimikosSolver_hh
#define solvers_StratimikosSolver_hh

#include <string>

#include <SPn/config.h>

#include "comm/P_Stream.hh"
#include "LinearSolver.hh"

#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_PreconditionerBase.hpp>

namespace profugus
{

//===========================================================================//
/*!
 * \class StratimikosSolver
 * \brief Use Trilinos' Stratimikos solver interface to solve a linear system.
 *
 * This class provides a general interface to Trilinos's Stratimikos linear
 * solver package.  The user provides a Thyra LinearOpBase and two Thyra
 * MultiVecBase objects, the following system is solved:
 * \f[
   \mathbf{A}X = B\:.
 * \f]
 */
/*!
 * \example solvers/test/tstStratimikos.cc
 *
 * Test of StratimikosSolver.
 */
//===========================================================================//

class StratimikosSolver :
    public LinearSolver<Epetra_MultiVector,Epetra_Operator>
{
  public:
    //@{
    //! Useful typedefs.
    typedef Teuchos::ParameterList                      ParameterList;
    typedef Teuchos::RCP<ParameterList>                 RCP_ParameterList;
    typedef Epetra_MultiVector                          MV;
    typedef Epetra_Operator                             OP;
    typedef Thyra::LinearOpBase<double>                 LOp;
    typedef Thyra::PreconditionerBase<double>           Prec;
    typedef Thyra::LinearOpWithSolveBase<double>        LOWS;
    typedef Thyra::LinearOpWithSolveFactoryBase<double> LOWS_Factory;
    typedef Teuchos::RCP<LOWS_Factory>                  RCP_LOWS_Factory;
    typedef Teuchos::RCP<const LOp>                     RCP_LOp;
    typedef Teuchos::RCP<LOWS>                          RCP_LOWS;
    typedef Teuchos::RCP<const Prec>                    RCP_Prec;
    //@}

  private:
    RCP_LOWS_Factory   d_factory;
    RCP_LOWS           d_solver;
    RCP_LOp            d_thyraA;
    RCP_LOp            d_prec;

  public:
    // Constructor.
    // Read Profugus database entries for solver parameters.
    explicit StratimikosSolver(RCP_ParameterList db);

    // Set Epetra Operator for linear system
    void set_operator(Teuchos::RCP<Epetra_Operator> A);

    // Set Epetra Operator for preconditioner
    void set_preconditioner(Teuchos::RCP<Epetra_Operator> P);

    // Solve a linear problem.
    void solve(Teuchos::RCP<Epetra_MultiVector>       ep_x,
               Teuchos::RCP<const Epetra_MultiVector> ep_b);

    // >>> ACCESSORS

    //! Set tolerance.
    double tolerance() const { return b_tolerance; }

    //! Maximum number of iterations.
    int max_itr() const { return b_max_iters; }
};

} // end namespace profugus

#endif // solvers_StratimikosSolver_hh

//---------------------------------------------------------------------------//
//              end of solvers/StratimikosSolver.hh
//---------------------------------------------------------------------------//
