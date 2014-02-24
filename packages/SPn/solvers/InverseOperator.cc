//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/InverseOperator.cc
 * \author Thomas Evans, Steven Hamilton
 * \date   Fri Feb 21 13:05:35 2014
 * \brief  InverseOperator member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "InverseOperator.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 */
InverseOperator::InverseOperator(RCP_ParameterList db)
{
    d_solver = LinearSolverBuilder::build_solver(db);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set epetra operator
 *
 * \param A Epetra_Operator
 */
void InverseOperator::set_operator( RCP_Operator A )
{
    Require( !A.is_null() );
    d_solver->set_operator(A);

    // Store A for some of the Epetra interface functions
    d_A = A;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set epetra operator for rhs
 *
 * \param B Epetra_Operator
 */
void InverseOperator::set_rhs_operator( RCP_Operator B )
{
    Require( !B.is_null() );
    d_B = B;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set epetra operator for preconditioner
 *
 * \param P Epetra_Operator
 */
void InverseOperator::set_preconditioner( RCP_Operator P )
{
    Require( !P.is_null() );
    Require( !d_solver.is_null() );
    d_solver->set_preconditioner(P);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform y=A^{-1}x or y=A^{-1}Bx
 *
 * \param x Input vector
 * \param y Output vector
 */
int InverseOperator::Apply(const MV &x,
                                 MV &y ) const
{
    Require( x.MyLength() == x.MyLength() );
    Require( d_A->OperatorDomainMap().NumMyElements()==x.MyLength() );

    if( !(d_B.is_null()) )
    {
        Require( d_B->OperatorDomainMap().NumMyElements()==x.MyLength() );
        MV z(x);
        d_B->Apply(x,z);
        d_solver->solve(Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(z));
    }
    else
    {
        d_solver->solve(Teuchos::rcpFromRef(y), Teuchos::rcpFromRef(x));
    }

    return 0;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of InverseOperator.cc
//---------------------------------------------------------------------------//
