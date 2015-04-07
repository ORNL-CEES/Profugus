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

#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_PreconditionerBase.hpp"

#include "AnasaziMultiVecTraits.hpp"

namespace profugus
{

//===========================================================================//
/*!
 * \class StratimikosSolver
 * \brief Use Trilinos' Stratimikos solver interface to solve a linear system.
 *
 * This class provides a general interface to Trilinos's Stratimikos linear
 * solver package.  The user provides a Epetra/Tpetra operators and
 * multivectors and this class wraps it into appropriate Thyra operators and
 * vectors, builds a solver using Stratimikos, and solves the system.
 */
/*!
 * \example solvers/test/tstStratimikos.cc
 *
 * Test of StratimikosSolver.
 */
//===========================================================================//

template <class T>
class StratimikosSolver : public LinearSolver<T>
{
  public:
    //@{
    //! Useful typedefs.
    typedef typename T::ST                              ST;
    typedef typename T::MV                              MV;
    typedef typename T::OP                              OP;
    typedef LinearSolver<T>                             Base;
    typedef Teuchos::ParameterList                      ParameterList;
    typedef Thyra::LinearOpBase<double>                 LOp;
    typedef Thyra::PreconditionerBase<double>           Prec;
    typedef Thyra::LinearOpWithSolveBase<double>        LOWS;
    typedef Thyra::LinearOpWithSolveFactoryBase<double> LOWS_Factory;
    typedef Anasazi::MultiVecTraits<double,MV>          MVT;
    //@}

  private:

    Teuchos::RCP<LOWS_Factory>  d_factory;
    Teuchos::RCP<LOWS>          d_solver;
    Teuchos::RCP<const LOp>     d_thyraA;
    Teuchos::RCP<const LOp>     d_prec;
    bool                        d_updated_operator;

    using Base::b_tolerance;
    using Base::b_verbosity;
    using Base::b_max_iters;
    using Base::b_num_iters;
    using Base::b_label;

  public:
    // Constructor.
    // Read Profugus database entries for solver parameters.
    explicit StratimikosSolver(Teuchos::RCP<ParameterList> db);

    // Set Operator for linear system
    void set_operator(Teuchos::RCP<OP> A);

    // Set Operator for preconditioner
    void set_preconditioner(Teuchos::RCP<OP> P);

    // Solve a linear problem.
    void solve(Teuchos::RCP<MV>       x,
               Teuchos::RCP<const MV> b);

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
