//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/MueLuPreconditioner.t.hh
 * \author Thomas M. Evans, Steven P. Hamilton
 * \date   Fri Feb 21 13:41:13 2014
 * \brief  MueLuPreconditioner template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_MueLuPreconditioner_t_hh
#define solvers_MueLuPreconditioner_t_hh

#include "MueLuPreconditioner.hh"

namespace profugus
{

//
// Implementation of MueLuPreconditionerBase
//

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 */
MueLuPreconditionerBase::MueLuPreconditionerBase()
{
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set up the MueLu problem
 *
 * \param pl List of problem parameters
 */
void MueLuPreconditionerBase::setup(
    Teuchos::RCP<Teuchos::ParameterList> pl )
{
    Require( pl != Teuchos::null );
    Require( d_matrix != Teuchos::null );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply MueLu preconditioner to vector
 *
 * \param x Input vector
 * \param y Output vector
 */
void MueLuPreconditionerBase::ApplyImpl(
        Teuchos::RCP<const Xpetra_MultiVector> x,
        Teuchos::RCP<Xpetra_MultiVector>       y ) const
{
}

} // end namespace profugus

#endif // solvers_MueLuPreconditioner_t_hh

//---------------------------------------------------------------------------//
//                 end of MueLuPreconditioner.t.hh
//---------------------------------------------------------------------------//
