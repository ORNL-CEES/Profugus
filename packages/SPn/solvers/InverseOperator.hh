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
#include "Teuchos_RCP.hpp"

#include "harness/DBC.hh"
#include "LinearSolver.hh"
#include "LinAlgTypedefs.hh"

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
template <class T>
class InverseOperatorBase
{
  public:

    typedef typename T::MV                        MV;
    typedef typename T::OP                        OP;
    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    typedef Anasazi::OperatorTraits<double,MV,OP> OPT;

    InverseOperatorBase(Teuchos::RCP<Teuchos::ParameterList> pl);
    virtual void set_operator( Teuchos::RCP<OP> A );
    virtual void set_rhs_operator( Teuchos::RCP<OP> B );
    void set_preconditioner( Teuchos::RCP<OP> P );

  protected:

    void ApplyImpl( const MV &x, MV &y ) const;

    Teuchos::RCP<OP> d_A, d_B, d_P;
    Teuchos::RCP<LinearSolver<T> > d_solver;
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
template <class T>
class InverseOperator
{
  public:

    InverseOperator();
};

// Implementation for Epetra_MultiVector/Operator
template <>
class InverseOperator<EpetraTypes>
    : public Epetra_Operator,
      public InverseOperatorBase<EpetraTypes>
{
  public:
    //@{
    //! Typedefs.
    typedef typename EpetraTypes::MV MV;
    typedef typename EpetraTypes::OP OP;
    typedef Teuchos::RCP<OP>         RCP_Operator;
    //@}

  private:

    typedef InverseOperatorBase<EpetraTypes> Base;
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
        REQUIRE( x.MyLength() == y.MyLength() );
        REQUIRE( d_A->OperatorDomainMap().NumMyElements()==x.MyLength() );

        if( !(d_B.is_null()) )
        {
            REQUIRE( d_B->OperatorDomainMap().NumMyElements()==x.MyLength() );
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
        REQUIRE(d_A != Teuchos::null);
        return d_A->Comm();
    }
    const Epetra_Map & OperatorDomainMap() const
    {
        if( d_B != Teuchos::null )
        {
            return d_B->OperatorDomainMap();
        }
        else
        {
            REQUIRE(d_A != Teuchos::null);
            return d_A->OperatorRangeMap();
        }
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        REQUIRE(d_A != Teuchos::null);
        return d_A->OperatorDomainMap();
    }

};

// Implementation for Tpetra::MultiVector/Operator
template <>
class InverseOperator<TpetraTypes>
    : public TpetraTypes::OP,
      public InverseOperatorBase<TpetraTypes>
{
  public:
    //@{
    //! Typedefs.
    typedef typename TpetraTypes::MV           MV;
    typedef typename TpetraTypes::OP           OP;
    typedef typename TpetraTypes::MAP          MAP;
    typedef Teuchos::RCP<OP>                   RCP_Operator;
    typedef Anasazi::MultiVecTraits<double,MV> MVT;
    //@}

  private:

    typedef InverseOperatorBase<TpetraTypes> Base;
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
        REQUIRE( x.getLocalLength() == y.getLocalLength() );
        REQUIRE( d_A->getDomainMap()->getNodeNumElements() ==
                 x.getLocalLength() );

        Teuchos::RCP<MV> z = MVT::CloneCopy(x);

        if( !(d_B.is_null()) )
        {
            REQUIRE( d_B->getDomainMap()->getNodeNumElements() ==
                     x.getLocalLength() );
        }

        ApplyImpl(x,y);

        if( beta == 0.0 )
        {
            y.scale(alpha);
        }
        else
        {
            y.update(beta,*z,alpha);
        }
    }

    // Required inherited interface.
    bool hasTranposeApply() const {return false;}
    Teuchos::RCP<const MAP> getDomainMap() const
    {
        REQUIRE(d_A != Teuchos::null);
        if( d_B != Teuchos::null )
        {
            return d_B->getDomainMap();
        }
        else
        {
            REQUIRE(d_A != Teuchos::null);
            return d_A->getRangeMap();
        }
    }
    Teuchos::RCP<const MAP> getRangeMap() const
    {
        REQUIRE(d_A != Teuchos::null);
        return d_A->getDomainMap();
    }

};

} // end namespace profugus

#endif // solvers_InverseOperator_hh

#include "InverseOperator.t.hh"

//---------------------------------------------------------------------------//
//                 end of InverseOperator.hh
//---------------------------------------------------------------------------//
