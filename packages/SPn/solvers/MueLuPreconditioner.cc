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
void MueLuPreconditionerBase::setup( Teuchos::RCP<Teuchos::ParameterList> pl )
{
    REQUIRE( pl != Teuchos::null );
    REQUIRE( d_matrix != Teuchos::null );

    Teuchos::RCP<MueLu::HierarchyManager<SCALAR,LO,GO,NODE> > mueLuFactory =
        Teuchos::rcp(
            new MueLu::ParameterListInterpreter<SCALAR,LO,GO,NODE>(*pl));

    d_hierarchy = mueLuFactory->CreateHierarchy();

    Teuchos::RCP<MueLu::Level> L = d_hierarchy->GetLevel(0);
    L->Set("A",d_matrix);

    mueLuFactory->SetupHierarchy(*d_hierarchy);

    d_hierarchy->IsPreconditioner(true);
}

} // end namespace profugus

#endif // solvers_MueLuPreconditioner_t_hh

//---------------------------------------------------------------------------//
//                 end of MueLuPreconditioner.t.hh
//---------------------------------------------------------------------------//
