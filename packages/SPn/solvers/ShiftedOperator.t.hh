//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/ShiftedOperator.t.hh
 * \author Thomas M. Evans, Steven P. Hamilton
 * \date   Fri Feb 21 13:41:13 2014
 * \brief  ShiftedOperator template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_solvers_ShiftedOperator_t_hh
#define SPn_solvers_ShiftedOperator_t_hh

#include "ShiftedOperator.hh"

namespace profugus
{

//
// Implementation of ShiftedOperatorBase
//

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 */
template <class T>
ShiftedOperatorBase<T>::ShiftedOperatorBase()
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
template <class T>
void ShiftedOperatorBase<T>::set_operator( Teuchos::RCP<OP> A )
{
    REQUIRE( !A.is_null() );
    d_A = A;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set epetra operator for rhs
 *
 * \param B Epetra_Operator
 */
template <class T>
void ShiftedOperatorBase<T>::set_rhs_operator( Teuchos::RCP<OP> B )
{
    REQUIRE( !B.is_null() );
    d_B = B;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform y=(A-lambda*I)x or y=(A-lambda*B)x
 *
 * \param x Input vector
 * \param y Output vector
 */
template <class T>
void ShiftedOperatorBase<T>::ApplyImpl(const MV &x,
                                                 MV &y ) const
{
    if( !(d_B.is_null()) )
    {
        Teuchos::RCP<MV> z = MVT::Clone(x,MVT::GetNumberVecs(x));
        OPT::Apply(*d_A,x,y);
        OPT::Apply(*d_B,x,*z);
        MVT::MvAddMv(-d_shift,*z,1.0,y,y);
    }
    else
    {
        OPT::Apply(*d_A,x,y);
        MVT::MvAddMv(-d_shift,x,1.0,y,y);
    }
}

} // end namespace profugus

#endif // SPn_solvers_ShiftedOperator_t_hh

//---------------------------------------------------------------------------//
//                 end of ShiftedOperator.t.hh
//---------------------------------------------------------------------------//
