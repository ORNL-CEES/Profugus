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

#include "MueLu_EasyParameterListInterpreter.hpp"

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

    // As of 9/4/14, MueLu throws an exception if constructed with default
    // parameters unless SuperLU is available.  We switch to a different
    // default coarse grid solver to avoid this.
    pl->get("coarse: type",std::string("RELAXATION"));

    // The template parameters here are a bit tricky...
    // We didn't want this class to be templated, but we need access
    // to SCALAR/LO/GO/NODE values for the HierarchyManager.  Right
    // now we're just pulling them from our Tpetra type which will
    // work as long as the Tpetra type is consistent with Epetra
    // (i.e. double/int/int).  If we change the Tpetra template parameters
    // this may cause problems with an Epetra instantiation, in which case
    // we would need to have this class be templated.
    Teuchos::RCP<MueLu::HierarchyManager<ST,LO,GO,NODE> > mueLuFactory =
        Teuchos::rcp(
            new MueLu::EasyParameterListInterpreter<ST,LO,GO,NODE>(*pl));

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
