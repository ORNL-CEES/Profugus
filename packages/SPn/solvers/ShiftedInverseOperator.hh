//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/ShiftedInverseOperator.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 13:38:42 2014
 * \brief  ShiftedInverseOperator class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_ShiftedInverseOperator_hh
#define solvers_ShiftedInverseOperator_hh

#include <Teuchos_RCP.hpp>

#include "InverseOperator.hh"
#include "ShiftedOperator.hh"

#include "LinAlgTypedefs.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class ShiftedInverseOperator
 * \brief Perform shift-and-invert operation.
 *
 * This class performs the operation \f$ y = (A - \lambda I)^{-1}x \f$,
 * or \f$ y = (A - \lambda B)^{-1} Bx \f$.
 * \brief
 */
/*!
 * \example solvers/test/tstShiftedInverseOperator.cc
 *
 * Test of ShiftedInverseOperator.
 */
//===========================================================================//

template <class T>
class ShiftedInverseOperator
{
    ShiftedInverseOperator();
};

template <>
class ShiftedInverseOperator<EpetraTypes>
    : public InverseOperator<EpetraTypes>
{
  private:
    typedef typename EpetraTypes::MV     MV;
    typedef typename EpetraTypes::OP     OP;
    typedef InverseOperator<EpetraTypes> Base;

    // >>> DATA
    Teuchos::RCP<ShiftedOperator<EpetraTypes> > d_operator;
    double d_shift;

  public:

    // Constructor.
    explicit ShiftedInverseOperator(Teuchos::RCP<Teuchos::ParameterList> pl)
        : Base(pl)
        , d_operator( Teuchos::rcp( new ShiftedOperator<EpetraTypes>() ) )
        , d_shift(0.0)
    {
    }

    // Set shift
    void set_shift( double shift )
    {
        d_operator->set_shift(shift);
    }

    // Set "A" operator
    void set_operator( Teuchos::RCP<OP> A )
    {
        d_operator->set_operator(A);
        Base::set_operator(d_operator);
    }

    // Set "B" operator
    void set_rhs_operator( Teuchos::RCP<OP> B )
    {
        d_operator->set_rhs_operator(B);
        Base::set_rhs_operator(B);
    }

    // Apply (solve linear system)
    int Apply(const MV &x, MV &y) const
    {
        REQUIRE( !d_operator.is_null() );
        REQUIRE( x.MyLength() == y.MyLength() );
        REQUIRE( d_operator->OperatorDomainMap().NumMyElements()==x.MyLength() );

        d_solver->solve(Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x));

        return 0;
    }
};

template <>
class ShiftedInverseOperator<TpetraTypes>
    : public InverseOperator<TpetraTypes>
{
  private:
    typedef typename TpetraTypes::MV           MV;
    typedef typename TpetraTypes::OP           OP;
    typedef InverseOperator<TpetraTypes>            Base;
    typedef Anasazi::MultiVecTraits<double,MV> MVT;

    // >>> DATA
    Teuchos::RCP<ShiftedOperator<TpetraTypes> > d_operator;
    double d_shift;

  public:

    // Constructor.
    explicit ShiftedInverseOperator(Teuchos::RCP<Teuchos::ParameterList> pl)
        : Base(pl)
        , d_operator( Teuchos::rcp( new ShiftedOperator<TpetraTypes>() ) )
        , d_shift(0.0)
    {}

    // Set shift
    void set_shift( double shift )
    {
        d_operator->set_shift(shift);
    }
    //
    // Set "A" operator
    void set_operator( Teuchos::RCP<OP> A )
    {
        d_operator->set_operator(A);
        Base::set_operator(d_operator);
    }

    // Set "B" operator
    void set_rhs_operator( Teuchos::RCP<OP> B )
    {
        d_operator->set_rhs_operator(B);
        Base::set_rhs_operator(B);
    }

    // Apply (solve linear system)
    void apply(const MV &x, MV &y,
               Teuchos::ETransp mode=Teuchos::NO_TRANS,
               double alpha=1.0, double beta=0.0) const
    {
        REQUIRE( !d_operator.is_null() );
        REQUIRE( x.getLocalLength() == y.getLocalLength() );
        REQUIRE( d_operator->getDomainMap()->getNodeNumElements() ==
                 x.getLocalLength() );

        Teuchos::RCP<MV> z = MVT::CloneCopy(x);

        d_solver->solve(Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x));

        if( beta == 0.0 )
        {
            y.scale(alpha);
        }
        else
        {
            y.update(beta,*z,alpha);
        }
    }
};

} // end namespace profugus

#endif // solvers_ShiftedInverseOperator_hh

//---------------------------------------------------------------------------//
//                 end of ShiftedInverseOperator.hh
//---------------------------------------------------------------------------//
