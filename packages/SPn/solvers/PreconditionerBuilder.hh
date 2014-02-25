//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/PreconditionerBuilder.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 09:28:52 2014
 * \brief  PreconditionerBuilder class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_PreconditionerBuilder_hh
#define solvers_PreconditionerBuilder_hh

#include "harness/DBC.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"

namespace profugus
{

//===========================================================================//
/*!
 * \class Preconditioner
 * \brief Persistence class for Epetra-based preconditioners.
 *
 * Storage for matrix-based preconditioners.  The reason for this class is
 * that Ifpack/ML objects are applied with the "ApplyInverse" function (as is
 * the Aztec convention) and our convention is to apply them with "Apply"
 * (this is the Anasazi/Belos approach).  The convention can be reversed by
 * wrapping the "raw" preconditioner inside an "Epetra_InvOperator."  However,
 * because Epetra_InvOperator only stores a raw pointer to the original
 * operator, a copy of that original must be maintained for the entire life of
 * the preconditioner.
 *
 */
//===========================================================================//

class Preconditioner : public Epetra_Operator
{
  public:
    // Constructor
    explicit Preconditioner( Teuchos::RCP<Epetra_Operator> raw_prec )
        : d_raw_prec(raw_prec)
    {
        Require( d_raw_prec != Teuchos::null );
    }

    // Apply
    virtual int Apply(const Epetra_MultiVector &x, Epetra_MultiVector &y) const
    {
        Require( d_raw_prec != Teuchos::null );
        return d_raw_prec->ApplyInverse(x,y);
    }

    // Required inherited interface.
    int SetUseTranspose(bool use_transpose)
    {
        return d_raw_prec->SetUseTranspose(use_transpose);
    }

    int ApplyInverse(const Epetra_MultiVector &x, Epetra_MultiVector &y) const
    {
        return d_raw_prec->Apply(x,y);
    }

    double NormInf() const { return d_raw_prec->NormInf(); }
    const char *  Label() const { return "Exnihilo Preconditioner"; }
    bool UseTranspose() const { return d_raw_prec->UseTranspose(); }
    bool HasNormInf() const { return d_raw_prec->HasNormInf(); }
    const Epetra_Comm & Comm() const
    {
        return d_raw_prec->Comm();
    }
    const Epetra_Map & OperatorDomainMap() const
    {
        return d_raw_prec->OperatorRangeMap();
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        return d_raw_prec->OperatorDomainMap();
    }

  private:

    Teuchos::RCP<Epetra_Operator> d_raw_prec;
};

//===========================================================================//
/*!
 * \class PreconditionerBuilder
 * \brief
 */
/*!
 * \example solvers/test/tstPreconditionerBuilder.cc
 *
 * Test of PreconditionerBuilder.
 */
//===========================================================================//

class PreconditionerBuilder
{
  public:

    typedef Teuchos::RCP<Teuchos::ParameterList> RCP_ParameterList;
    typedef Epetra_MultiVector                   MV;
    typedef Epetra_Operator                      OP;

    static Teuchos::RCP<Epetra_Operator> build_preconditioner(
        Teuchos::RCP<Epetra_Operator> op, RCP_ParameterList db );
};

} // end namespace profugus

#endif // solvers_PreconditionerBuilder_hh

//---------------------------------------------------------------------------//
//                 end of PreconditionerBuilder.hh
//---------------------------------------------------------------------------//
