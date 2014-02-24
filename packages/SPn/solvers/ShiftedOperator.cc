//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/ShiftedOperator.cc
 * \author Thomas M. Evans, Steven P. Hamilton
 * \date   Fri Feb 21 13:41:13 2014
 * \brief  ShiftedOperator member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ShiftedOperator.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 */
ShiftedOperator::ShiftedOperator()
    : d_shift(0.0)
{
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set epetra operator
 *
 * \param A Epetra_Operator
 */
void ShiftedOperator::set_operator( RCP_Operator A )
{
    Require( !A.is_null() );
    d_A = A;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set epetra operator for rhs
 *
 * \param B Epetra_Operator
 */
void ShiftedOperator::set_rhs_operator( RCP_Operator B )
{
    Require( !B.is_null() );
    d_B = B;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform y=(A-lambda*I)x or y=(A-lambda*B)x
 *
 * \param x Input vector
 * \param y Output vector
 */
int ShiftedOperator::Apply(const MV &x,
                                 MV &y ) const
{
    Require( x.MyLength() == y.MyLength() );
    Require( d_A->OperatorDomainMap().NumMyElements()==x.MyLength() );

    if( !(d_B.is_null()) )
    {
        Require( d_B->OperatorDomainMap().NumMyElements()==x.MyLength() );
        MV z(x);
        d_A->Apply(x,y);
        d_B->Apply(x,z);
        y.Update(-d_shift,z,1.0);
    }
    else
    {
        d_A->Apply(x,y);
        y.Update(-d_shift,x,1.0);
    }

    return 0;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of ShiftedOperator.cc
//---------------------------------------------------------------------------//
