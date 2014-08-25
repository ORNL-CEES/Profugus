//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/InverseOperator.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 13:41:13 2014
 * \brief  InverseOperator class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_InverseOperator_hh
#define solvers_InverseOperator_hh

#include <string>

#include <SPn/config.h>

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"

#include "harness/DBC.hh"
#include "LinearSolver.hh"

namespace profugus
{

//===========================================================================//
/*
 *! \class InverseOperatorBase
 *  \brief Implementation details of InverseOperator.
 *
 *  This class performs the operation y = A^{-1}x or y = A^{-1}Bx.
 *  This class should not be constructed directly, but is rather an
 *  implementation detail of the InverseOperator class, which will
 *  additionally implement either the Epetra or Tpetra operator interface.
 */
//===========================================================================//
template <class MV, class OP>
class InverseOperatorBase
{
  public:

    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    typedef Anasazi::OperatorTraits<double,MV,OP> OPT;

    InverseOperatorBase(Teuchos::RCP<Teuchos::ParameterList> pl);
    void set_operator( Teuchos::RCP<OP> A );
    void set_rhs_operator( Teuchos::RCP<OP> B );
    void set_preconditioner( Teuchos::RCP<OP> P );

  protected:

    void ApplyImpl( const MV &x, MV &y ) const;

    Teuchos::RCP<OP> d_A, d_B, d_P;
    Teuchos::RCP<LinearSolver<MV,OP> > d_solver;
};

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
//

// Dummy implementation for MV/OP combos lacking a specialization.
// Attempt to instantiate will cause a compile error.
template <class MV, class OP>
class InverseOperator
{
  public:

    InverseOperator();
};

// Implementation for Epetra_MultiVector/Operator
template <>
class InverseOperator<Epetra_MultiVector,Epetra_Operator>
    : public Epetra_Operator,
      public InverseOperatorBase<Epetra_MultiVector,Epetra_Operator>
{
  public:
    //@{
    //! Typedefs.
    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator    OP;
    typedef Teuchos::RCP<OP>   RCP_Operator;
    //@}

  private:

    typedef InverseOperatorBase<MV,OP> Base;
    using Base::d_A;
    using Base::d_B;

  public:

    // Constructor.
    explicit InverseOperator(Teuchos::RCP<Teuchos::ParameterList> pl)
        : Base(pl)
    {}

    // Apply (solve linear system)
    int Apply(const MV &x, MV &y) const
    {
        Require( x.MyLength() == y.MyLength() );
        Require( d_A->OperatorDomainMap().NumMyElements()==x.MyLength() );

        if( !(d_B.is_null()) )
        {
            Require( d_B->OperatorDomainMap().NumMyElements()==x.MyLength() );
        }

        ApplyImpl(x,y);

        return 0;
    }

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
    const Epetra_Map & OperatorDomainMap() const
    {
        Require(d_A != Teuchos::null);
        if( d_B != Teuchos::null )
        {
            return d_B->OperatorDomainMap();
        }
        else
        {
            Require(d_A != Teuchos::null);
            return d_A->OperatorRangeMap();
        }
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        Require(d_A != Teuchos::null);
        return d_A->OperatorDomainMap();
    }

};

// Implementation for Tpetra::MultiVector/Operator
template <>
class InverseOperator<
        Tpetra::MultiVector<double,int,int,KokkosClassic::SerialNode>,
        Tpetra::Operator<double,int,int,KokkosClassic::SerialNode> >
    : public Tpetra::Operator<double,int,int,KokkosClassic::SerialNode>,
      public InverseOperatorBase<
                Tpetra::MultiVector<double,int,int,KokkosClassic::SerialNode>,
                Tpetra::Operator<double,int,int,KokkosClassic::SerialNode> >
{
  public:
    //@{
    //! Typedefs.
    typedef Tpetra::MultiVector<double,int,int,KokkosClassic::SerialNode> MV;
    typedef Tpetra::Operator<double,int,int,KokkosClassic::SerialNode>    OP;
    typedef Tpetra::Map<int,int,KokkosClassic::SerialNode>                Map;
    typedef Teuchos::RCP<OP>   RCP_Operator;
    //@}

  private:

    typedef InverseOperatorBase<MV,OP> Base;
    using Base::d_A;
    using Base::d_B;

  public:

    // Constructor.
    explicit InverseOperator(Teuchos::RCP<Teuchos::ParameterList> pl)
        : Base(pl)
    {};

    // Apply (solve linear system)
    void apply(const MV &x, MV &y,
               Teuchos::ETransp mode=Teuchos::NO_TRANS,
               double alpha=1.0, double beta=0.0) const
    {
        Require( alpha = 1.0 );
        Require( beta = 0.0 );
        Require( x.getLocalLength() == y.getLocalLength() );
        Require( d_A->getDomainMap()->getNodeNumElements() ==
                 x.getLocalLength() );

        if( !(d_B.is_null()) )
        {
            Require( d_B->getDomainMap()->getNodeNumElements() ==
                     x.getLocalLength() );
        }

        ApplyImpl(x,y);
    }

    // Required inherited interface.
    bool hasTranposeApply() const {return false;}
    Teuchos::RCP<const Map> getDomainMap() const
    {
        Require(d_A != Teuchos::null);
        if( d_B != Teuchos::null )
        {
            return d_B->getDomainMap();
        }
        else
        {
            Require(d_A != Teuchos::null);
            return d_A->getRangeMap();
        }
    }
    Teuchos::RCP<const Map> getRangeMap() const
    {
        Require(d_A != Teuchos::null);
        return d_A->getDomainMap();
    }

};

} // end namespace profugus

#endif // solvers_InverseOperator_hh

#include "InverseOperator.t.hh"

//---------------------------------------------------------------------------//
//                 end of InverseOperator.hh
//---------------------------------------------------------------------------//
