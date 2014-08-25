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

namespace profugus
{

//===========================================================================//
/*
 *! \class ShiftedOperatorBase
 *  \brief Implementation details of ShiftedOperator.
 *
 *  This class performs the operation y = (A - \lambda B)*x
 *  This class should not be constructed directly, but is rather an
 *  implementation detail of the ShiftedOperator class, which will
 *  additionally implement either the Epetra or Tpetra operator interface.
 */
//===========================================================================//
template <class MV, class OP>
class ShiftedOperatorBase
{
  public:

    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    typedef Anasazi::OperatorTraits<double,MV,OP> OPT;

    ShiftedOperatorBase();
    void set_operator( Teuchos::RCP<OP> A );
    void set_rhs_operator( Teuchos::RCP<OP> B );
    void set_shift(double shift){ d_shift = shift; }

  protected:

    void ApplyImpl( const MV &x, MV &y ) const;

    Teuchos::RCP<OP> d_A, d_B;
    double d_shift;
};

//===========================================================================//
/*!
 * \class ShiftedOperator
 * \brief Perform shifted apply (A - \lambda B)*x
 *
 * This class takes two linear operators and allows the application of
 * the shifted operator (A - \lambda B)*x as a single operator as an
 * Epetra/Tpetra operator.  If B is not provided, the identity will be
 * used, producing the operator (A - \lambda I)*x.
 */
/*!
 * \example solvers/test/tstShiftedOperator.cc
 *
 * Test of ShiftedOperator.
 */
//===========================================================================//

// Dummy implementation for MV/OP combos lacking a specialization.
// Attempt to instantiate will cause a compile error.
template <class MV, class OP>
class ShiftedOperator
{
  public:

    ShiftedOperator();
};

// Implementation for Epetra_MultiVector/Operator
template <>
class ShiftedOperator<Epetra_MultiVector,Epetra_Operator>
    : public Epetra_Operator,
      public ShiftedOperatorBase<Epetra_MultiVector,Epetra_Operator>
{
  public:
    //@{
    //! Typedefs.
    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator    OP;
    typedef Teuchos::RCP<OP>   RCP_Operator;
    //@}

  private:

    typedef ShiftedOperatorBase<MV,OP> Base;
    using Base::d_A;
    using Base::d_B;
    using Base::d_shift;

  public:

    // Constructor.
    explicit ShiftedOperator(){};

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

// Implementation for Tpetra::MultiVector/Operator
template <>
class ShiftedOperator<
        Tpetra::MultiVector<double,int,int,KokkosClassic::SerialNode>,
        Tpetra::Operator<double,int,int,KokkosClassic::SerialNode> >
    : public Tpetra::Operator<double,int,int,KokkosClassic::SerialNode>,
      public ShiftedOperatorBase<
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

    typedef ShiftedOperatorBase<MV,OP> Base;
    using Base::d_A;
    using Base::d_B;
    using Base::d_shift;

  public:

    // Constructor.
    explicit ShiftedOperator(){};

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

#endif // solvers_ShiftedOperator_hh

#include "ShiftedOperator.t.hh"

//---------------------------------------------------------------------------//
//                 end of ShiftedOperator.hh
//---------------------------------------------------------------------------//
