//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/ShiftedInverseOperator.cc
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 13:38:42 2014
 * \brief  ShiftedInverseOperator member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "ShiftedInverseOperator.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 */
ShiftedInverseOperator::ShiftedInverseOperator(RCP_ParameterList db)
    : Base(db)
    , d_operator( Teuchos::rcp(new ShiftedOperator) )
    , d_shift(0.0)
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
void ShiftedInverseOperator::set_operator( RCP_Operator A )
{
    Require( !A.is_null() );
    d_A = A;
    d_operator->set_operator(A);

    // Even though ShiftedOperator object exists at construction time,
    //  we can't give it to the solver until we have A because its
    //  Epetra interface relies on A
    d_solver->set_operator(d_operator);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set epetra operator for rhs
 *
 * \param B Epetra_Operator
 */
void ShiftedInverseOperator::set_rhs_operator( RCP_Operator B )
{
    Require( !B.is_null() );
    d_B = B;
    d_operator->set_rhs_operator(B);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform y=(A-lambda*I)^{-1}x or y=(A-lambda*B)^{-1}x
 *
 * \param x Input vector
 * \param y Output vector
 */
int ShiftedInverseOperator::Apply(const MV &x,
                                  MV       &y ) const
{
    Require( !d_operator.is_null() );
    Require( x.MyLength() == x.MyLength() );
    Require( d_operator->OperatorDomainMap().NumMyElements()==x.MyLength() );

    d_solver->solve(Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x));

    return 0;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of ShiftedInverseOperator.cc
//---------------------------------------------------------------------------//
