//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/InverseOperator.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 13:05:35 2014
 * \brief  InverseOperator class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_InverseOperator_hh
#define solvers_InverseOperator_hh

#include <string>

#include <SPn/config.h>

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"

#include "harness/DBC.hh"
#include "LinearSolverBuilder.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class InverseOperator
 * \brief Wrap LinearSolver into an Epetra_Operator interface
 *
 * This class wraps the LinearSolver class, which solves a linear
 * system \f$ Ax = b \f$, into an Epetra_Operator such that
 * the result of this linear solve can be obtained as Apply(b,x).
 * This operator can be roughly viewed as applying the inverse of A to b.
 * The primary purpose for this capability is to use a solver as a
 * preconditioner (for instance with a Davidson eigensolver).
 *
 * Additionally, if a second operator (B) is provided, it will be applied
 * to the provided vector before calling the solver.  This allows the
 * operation \f$ y = A^{-1}Bx \f$ to be performed as a single operator
 * this is useful for Arnoldi or Power Iteration style eigensovlers that
 * require converting a generalized eigenvalue problem to a standard
 * eigenvalue problem.
 */
/*!
 * \example solvers/test/tstInverseOperator.cc
 *
 * Test of InverseOperator.
 */
//===========================================================================//

class InverseOperator: public Epetra_Operator
{
  public:
    //@{
    //! Typedefs.
    typedef Epetra_MultiVector                            MV;
    typedef Epetra_Operator                               OP;
    typedef Teuchos::RCP<OP>                              RCP_Operator;
    typedef LinearSolverBuilder<MV,OP>::RCP_ParameterList RCP_ParameterList;
    typedef LinearSolverBuilder<MV,OP>::RCP_LinearSolver  RCP_LinearSolver;
    //@}

  protected:
    // >>> SHARED DATA

    RCP_Operator     d_A;
    RCP_Operator     d_B;
    RCP_LinearSolver d_solver;

  public:
    // Constructor.
    // Read database entries for solver parameters.
    explicit InverseOperator(RCP_ParameterList db);

    // Set Epetra Operator
    virtual void set_operator(RCP_Operator A);
    virtual void set_preconditioner(RCP_Operator P);
    virtual void set_rhs_operator(RCP_Operator B);

    // Apply (solve linear system)
    virtual int Apply(const MV &x, MV &y) const;

    // Required inherited interface.
    int SetUseTranspose(bool useTranspose){return -1;}
    int ApplyInverse(const MV &x, MV &b) const
    {
        return -1;
    }
    double NormInf() const { return 0.0; }
    const char *  Label() const { return "InverseOperator"; }
    bool UseTranspose() const { return false; }
    bool HasNormInf() const { return false; }
    const Epetra_Comm & Comm() const
    {
        Require(d_A != Teuchos::null);
        return d_A->Comm();
    }
    // Because this is an inverse map of A, the notion of
    //  domain and range are reversed
    const Epetra_Map & OperatorDomainMap() const
    {
        Require(d_A != Teuchos::null);
        return d_A->OperatorRangeMap();
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        Require(d_A != Teuchos::null);
        return d_A->OperatorDomainMap();
    }
};

} // end namespace profugus

#endif // solvers_InverseOperator_hh

//---------------------------------------------------------------------------//
//                 end of InverseOperator.hh
//---------------------------------------------------------------------------//
