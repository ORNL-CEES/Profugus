//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/InverseOperator.t.hh
 * \author Thomas Evans, Steven Hamilton
 * \date   Fri Feb 21 13:05:35 2014
 * \brief  InverseOperator template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_InverseOperator_t_hh
#define solvers_InverseOperator_t_hh

#include "InverseOperator.hh"
#include "LinearSolverBuilder.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 */
template <LinAlgType T>
InverseOperatorBase<T>::InverseOperatorBase(
    Teuchos::RCP<Teuchos::ParameterList> pl )
{
    d_solver = LinearSolverBuilder<T>::build_solver(pl);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set epetra operator
 *
 * \param A Epetra_Operator
 */
template <LinAlgType T>
void InverseOperatorBase<T>::set_operator( Teuchos::RCP<OP> A )
{
    REQUIRE( !A.is_null() );
    d_solver->set_operator(A);

    // Store A for some of the operator interface functions
    d_A = A;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set epetra operator for rhs
 *
 * \param B Epetra_Operator
 */
template <LinAlgType T>
void InverseOperatorBase<T>::set_rhs_operator( Teuchos::RCP<OP> B )
{
    REQUIRE( !B.is_null() );
    d_B = B;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set epetra operator for preconditioner
 *
 * \param P Epetra_Operator
 */
template <LinAlgType T>
void InverseOperatorBase<T>::set_preconditioner( Teuchos::RCP<OP> P )
{
    REQUIRE( !P.is_null() );
    REQUIRE( !d_solver.is_null() );
    d_solver->set_preconditioner(P);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform y=A^{-1}x or y=A^{-1}Bx
 *
 * \param x Input vector
 * \param y Output vector
 */
template <LinAlgType T>
void InverseOperatorBase<T>::ApplyImpl(const MV &x, MV &y ) const
{
    if( !(d_B.is_null()) )
    {
        Teuchos::RCP<MV> z = MVT::Clone(x,MVT::GetNumberVecs(x));
        OPT::Apply(*d_B,x,*z);
        d_solver->solve(Teuchos::rcpFromRef(y), z);
    }
    else
    {
        d_solver->solve(Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x));
    }
}

} // end namespace profugus

#endif // solvers_InverseOperator_t_hh

//---------------------------------------------------------------------------//
//                 end of InverseOperator.t.hh
//---------------------------------------------------------------------------//
