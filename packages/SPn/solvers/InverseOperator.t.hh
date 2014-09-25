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
template <class T>
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
template <class T>
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
template <class T>
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
template <class T>
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
template <class T>
void InverseOperatorBase<T>::ApplyImpl(const MV &x, MV &y ) const
{
    // The data in the output vector y will be used as an initial guess
    //  to the linear solver, which is a good idea if we're using this
    //  inside of power iteration, in which case y will be a very good
    //  approximation to the solution when the eigensolver is nearly
    //  converged.  However, if this appears inside Arnoldi, y may
    //  have uninitialized data in it (thanks Anasazi), which could
    //  cause problems for the linear solver.  So we try to identify
    //  "bad" data on input and initialize y to zero in this case.
    //  We can't just look for inf or nan data, values that are "large"
    //  (i.e. 10^80) will cause problems for the linear solver.  We set
    //  an arbitrary threshold and if the vector norm exceeds this value
    //  we'll reset it to zero.  This is a bit of a hack but I don't have
    //  a better alternative to support our use cases.
    double big_value = 1e10;
    std::vector<double> nrm(1);
    MVT::MvNorm(y,nrm);
    if( Teuchos::ScalarTraits<typename T::ST>::isnaninf(nrm[0]) ||
        nrm[0] > big_value )
    {
        MVT::MvInit(y,0.0);
    }

    // If we have a B, apply it before solving
    if( !(d_B.is_null()) )
    {
        Teuchos::RCP<MV> z = MVT::Clone(x,MVT::GetNumberVecs(x));
        OPT::Apply(*d_B,x,*z);
        d_solver->solve(Teuchos::rcpFromRef(y), z);
    }
    // Otherwise just solve
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
