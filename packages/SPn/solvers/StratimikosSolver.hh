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

template <LinAlgType T>
class StratimikosSolver : public LinearSolver<T>
{
  public:
    //@{
    //! Useful typedefs.
    typedef typename LinAlgTypedefs<T>::ST              ST;
    typedef typename LinAlgTypedefs<T>::LO              LO;
    typedef typename LinAlgTypedefs<T>::GO              GO;
    typedef typename LinAlgTypedefs<T>::NODE            NODE;
    typedef typename LinAlgTypedefs<T>::MV              MV;
    typedef typename LinAlgTypedefs<T>::OP              OP;
    typedef LinearSolver<T>                             Base;
    typedef Teuchos::ParameterList                      ParameterList;
    typedef Teuchos::RCP<ParameterList>                 RCP_ParameterList;
    typedef Thyra::LinearOpBase<double>                 LOp;
    typedef Thyra::PreconditionerBase<double>           Prec;
    typedef Thyra::LinearOpWithSolveBase<double>        LOWS;
    typedef Thyra::LinearOpWithSolveFactoryBase<double> LOWS_Factory;
    typedef Teuchos::RCP<LOWS_Factory>                  RCP_LOWS_Factory;
    typedef Teuchos::RCP<const LOp>                     RCP_LOp;
    typedef Teuchos::RCP<LOWS>                          RCP_LOWS;
    typedef Teuchos::RCP<const Prec>                    RCP_Prec;
    typedef Anasazi::MultiVecTraits<double,MV>          MVT;
    //@}

  private:
    RCP_LOWS_Factory   d_factory;
    RCP_LOWS           d_solver;
    RCP_LOp            d_thyraA;
    RCP_LOp            d_prec;
    bool               d_updated_operator;

    using Base::b_tolerance;
    using Base::b_verbosity;
    using Base::b_max_iters;
    using Base::b_num_iters;
    using Base::b_label;

  public:
    // Constructor.
    // Read Profugus database entries for solver parameters.
    explicit StratimikosSolver(RCP_ParameterList db);

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

  private:

    Teuchos::RCP<Thyra::MultiVectorBase<double> > buildThyraMV(
            Teuchos::RCP<MV> x,
            Teuchos::RCP<const Thyra::VectorSpaceBase<double> > space) const;

    Teuchos::RCP<const Thyra::MultiVectorBase<double> > buildThyraConstMV(
            Teuchos::RCP<const MV> x,
            Teuchos::RCP<const Thyra::VectorSpaceBase<double> > space) const;
};

} // end namespace profugus

#endif // solvers_StratimikosSolver_hh

//---------------------------------------------------------------------------//
//              end of solvers/StratimikosSolver.hh
//---------------------------------------------------------------------------//
