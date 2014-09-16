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
#include "Teuchos_RCP.hpp"

#include "LinAlgTypedefs.hh"

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
template <class T>
class ShiftedOperatorBase
{
  public:

    typedef typename T::MV                        MV;
    typedef typename T::OP                        OP;
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
template <class T>
class ShiftedOperator
{
  public:

    ShiftedOperator();
};

// Implementation for Epetra_MultiVector/Operator
template <>
class ShiftedOperator<EpetraTypes>
    : public EpetraTypes::OP,
      public ShiftedOperatorBase<EpetraTypes>
{
  public:
    //@{
    //! Typedefs.
    typedef typename EpetraTypes::MV MV;
    typedef typename EpetraTypes::OP OP;
    typedef Teuchos::RCP<OP>         RCP_Operator;
    //@}

  private:

    typedef ShiftedOperatorBase<EpetraTypes> Base;
    using Base::d_A;
    using Base::d_B;
    using Base::d_shift;

  public:

    // Constructor.
    explicit ShiftedOperator(){};

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
    const char *  Label() const { return "ShiftedOperator"; }
    bool UseTranspose() const { return false; }
    bool HasNormInf() const { return false; }
    const Epetra_Comm & Comm() const
    {
        REQUIRE(d_A != Teuchos::null);
        return d_A->Comm();
    }
    const Epetra_Map & OperatorDomainMap() const
    {
        REQUIRE(d_A != Teuchos::null);
        return d_A->OperatorDomainMap();
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        REQUIRE(d_A != Teuchos::null);
        return d_A->OperatorRangeMap();
    }

};

// Implementation for Tpetra::MultiVector/Operator
template <>
class ShiftedOperator<TpetraTypes>
    : public TpetraTypes::OP,
      public ShiftedOperatorBase<TpetraTypes>
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

    typedef ShiftedOperatorBase<TpetraTypes> Base;
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
        REQUIRE( x.getLocalLength() == y.getLocalLength() );
        REQUIRE( d_A->getDomainMap()->getNodeNumElements() ==
                 x.getLocalLength() );

        if( !(d_B.is_null()) )
        {
            REQUIRE( d_B->getDomainMap()->getNodeNumElements() ==
                     x.getLocalLength() );
        }

        Teuchos::RCP<MV> z = MVT::CloneCopy(x);

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

#endif // solvers_ShiftedOperator_hh

#include "ShiftedOperator.t.hh"

//---------------------------------------------------------------------------//
//                 end of ShiftedOperator.hh
//---------------------------------------------------------------------------//
