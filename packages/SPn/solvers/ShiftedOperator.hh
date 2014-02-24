//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/ShiftedOperator.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 13:41:13 2014
 * \brief  ShiftedOperator class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_ShiftedOperator_hh
#define solvers_ShiftedOperator_hh

#include <string>

#include <SPn/config.h>

#include <Epetra_Operator.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

#include "harness/DBC.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class ShiftedOperator
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
 * \example solvers/test/tstShiftedOperator.cc
 *
 * Test of ShiftedOperator.
 */
//===========================================================================//

class ShiftedOperator : public Epetra_Operator
{
  public:
    //@{
    //! Typedefs.
    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator    OP;
    typedef Teuchos::RCP<OP>   RCP_Operator;
    //@}

  private:
    // >>> DATA

    RCP_Operator d_A;
    RCP_Operator d_B;
    double d_shift;

  public:

    // Constructor.
    // Read Denovo database entries for solver parameters.
    explicit ShiftedOperator();

    // Set Epetra Operator
    void set_operator(RCP_Operator A);
    void set_rhs_operator(RCP_Operator B);
    void set_shift(double shift)
    {
        d_shift = shift;
    }

    // Apply (solve linear system)
    int Apply(const MV &x, MV &y) const;

    // Required inherited interface.
    int SetUseTranspose(bool useTranspose){return -1;}
    int ApplyInverse(const MV &x, MV &b) const
    {
        return -1;
    }
    double NormInf() const { return 0.0; }
    const char *  Label() const { return "ShiftedOperator"; }
    bool UseTranspose() const { return false; }
    bool HasNormInf() const { return false; }
    const Epetra_Comm & Comm() const
    {
        Require(d_A != Teuchos::null);
        return d_A->Comm();
    }
    const Epetra_Map & OperatorDomainMap() const
    {
        Require(d_A != Teuchos::null);
        return d_A->OperatorDomainMap();
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        Require(d_A != Teuchos::null);
        return d_A->OperatorRangeMap();
    }
};

} // end namespace profugus

#endif // solvers_ShiftedOperator_hh

//---------------------------------------------------------------------------//
//                 end of ShiftedOperator.hh
//---------------------------------------------------------------------------//
