//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/PreconditionerBuilder.t.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 09:28:52 2014
 * \brief  PreconditionerBuilder template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_PreconditionerBuilder_t_hh
#define solvers_PreconditionerBuilder_t_hh

#include <string>

#include <SPn/config.h>

#include "Epetra_RowMatrix.h"
#include "Epetra_InvOperator.h"
#include "Ifpack.h"

// ML has to be optional for Windows compatibility
#ifdef USE_ML
#include "ml_MultiLevelPreconditioner.h"
#endif

#include "utils/String_Functions.hh"
#include "PreconditionerBuilder.hh"
#include "MueLuPreconditioner.hh"

#include "Tpetra_Operator.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Ifpack2_Factory_decl.hpp"
#include "Ifpack2_Factory_def.hpp"

#include "TpetraTypedefs.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Build a preconditioner.
 *
 * This function creates a preconditioner from an Epetra_Operator.  It
 * is necessary that the underlying type of the operator is an
 * Epetra_RowMatrix.  Currently available preconditioner types (specified
 * via the "Preconditioner" db entry) are "Ifpack", "ML", and "None"
 */
template <>
Teuchos::RCP<Epetra_Operator>
PreconditionerBuilder<Epetra_Operator>::build_preconditioner(
    Teuchos::RCP<Epetra_Operator> op,
    RCP_ParameterList             db )
{
    using std::string;

    // Default to Ifpack
    string prec_type = to_lower(db->get("Preconditioner", string("ifpack")));
    VALIDATE(prec_type == "ifpack" || prec_type=="ml" ||
             prec_type=="none",
             "Preconditioner must be 'Ifpack', 'ML', or 'None'.");

    // Get the underlying matrix, must be Epetra_RowMatrix
    Teuchos::RCP<Epetra_RowMatrix> rowmat =
        Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(op);
    INSIST( rowmat != Teuchos::null,
            "Preconditioner construction requires Epetra_RowMatrix." );

    Teuchos::RCP<Epetra_Operator> prec;
    if( prec_type == "ifpack" )
    {
        // Create Ifpack preconditioner
        Ifpack ifpack_factory;
        Teuchos::RCP<Ifpack_Preconditioner> ifpack_prec;

        // Default to 0-overlap ILU
        string ifpack_type = db->get("Ifpack Type", string("ILU"));
        int overlap        = db->get("Ifpack Overlap", 0);
        ifpack_prec = Teuchos::rcp( ifpack_factory.Create(
            ifpack_type, rowmat.get(), overlap ) );

        // Convert "Ifpack Params" database to Teuchos::ParameterList
        RCP_ParameterList ifpack_db = Teuchos::sublist(db, "Ifpack Params");
        ifpack_prec->SetParameters(*ifpack_db);

        // Process preconditioner
        int err;
        err = ifpack_prec->Initialize();
        ENSURE( err == 0 );
        err = ifpack_prec->Compute();
        ENSURE( err == 0 );

        // Wrap raw preconditioner into an Epetra_InvOperator
        // to reverse the sense of Apply and ApplyInverse
        // Because Epetra_InvOperator stores a raw pointer rather than
        // an RCP, we need to keep the "raw" preconditioner alive.
        // We accomplish this by attaching the rcp to the raw
        // preconditioner as extra data on the inverted operator.
        prec = Teuchos::RCP<Epetra_Operator>(
            new Epetra_InvOperator(ifpack_prec.getRawPtr()) );
        Teuchos::set_extra_data(ifpack_prec,"ifpack_raw_pointer",Teuchos::inOutArg(prec));
        ENSURE( prec != Teuchos::null );
    }
    else if( prec_type == "ml" )
    {
#ifdef USE_ML
        Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ml_prec;

        // Create default db
        RCP_ParameterList ml_db = Teuchos::sublist(db, "ML Params");

        // Choose default settings
        string default_type = db->get("ML Default Type", string("DD"));

        // Set default profile, but don't override existing entries
        std::vector<int> az_options(AZ_OPTIONS_SIZE);
        std::vector<double> az_params(AZ_PARAMS_SIZE);
        bool override = false;
        ML_Epetra::SetDefaults(default_type, *ml_db, &az_options[0],
                               &az_params[0],override);

        ml_prec = Teuchos::rcp( new ML_Epetra::MultiLevelPreconditioner(
                                    *rowmat, *ml_db ) );

        // Wrap raw preconditioner into an Epetra_InvOperator
        // to reverse the sense of Apply and ApplyInverse
        // Because Epetra_InvOperator stores a raw pointer rather than
        // an RCP, we need to keep the "raw" preconditioner alive.
        // We accomplish this by attaching the rcp to the raw
        // preconditioner as extra data on the inverted operator.
        prec = Teuchos::RCP<Epetra_Operator>(
            new Epetra_InvOperator(ml_prec.getRawPtr()) );
        Teuchos::set_extra_data(ml_prec,"ml_raw_pointer",Teuchos::inOutArg(prec));

        ENSURE( prec != Teuchos::null );
#else
        VALIDATE(false,"ML not enabled in this build.");
#endif
    }

    return prec;
}

template <>
Teuchos::RCP<Tpetra_Operator>
PreconditionerBuilder<Tpetra_Operator>::build_preconditioner(
    Teuchos::RCP<Tpetra_Operator> op,
    RCP_ParameterList             db )
{
    string prec_type = to_lower(db->get("Preconditioner", string("none")));
    Teuchos::RCP<Tpetra_Operator> prec;
    if( prec_type == "ifpack2" )
    {
        // Dynamic cast to CrsMatrix
        Teuchos::RCP<Tpetra_CrsMatrix> row_mat =
            Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>( op );
        REQUIRE( row_mat != Teuchos::null );

        std::string ifpack2_type = db->get("Ifpack2_Type","ILUT");
        int overlap = db->get("Ifpack2_Overlap",0);

        Ifpack2::Factory factory;
        Teuchos::RCP<Teuchos::ParameterList> ifpack2_pl =
            Teuchos::sublist(db, "Ifpack2 Params");
        Teuchos::RCP<Ifpack2::Preconditioner<SCALAR,LO,GO,NODE> >
            ifpack_prec = factory.create(ifpack2_type,row_mat.getConst(),overlap);
        ifpack_prec->setParameters(*ifpack2_pl);
        ifpack_prec->initialize();
        ifpack_prec->compute();
        prec = ifpack_prec;
    }
    else if( prec_type == "muelu" )
    {
        // Dynamic cast to CrsMatrix
        Teuchos::RCP<Tpetra_CrsMatrix> row_mat =
            Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>( op );
        REQUIRE( row_mat != Teuchos::null );

        Teuchos::RCP<Teuchos::ParameterList> muelu_pl =
            Teuchos::sublist(db, "MueLu Params");

        // Wrap Tpetra objects as Xpetra
        prec = Teuchos::rcp(new MueLuPreconditioner<Tpetra_MultiVector,
                                                    Tpetra_Operator>(row_mat,
                                                                     muelu_pl));

    }
    else if( prec_type != "none" )
    {
        std::stringstream ss;
        ss << "Preconditioner " << prec_type << " not implemented" << std::endl;
        VALIDATE(false,ss.str());
    }
    return prec;
}

} // end namespace profugus

#endif // solvers_PreconditionerBuilder_t_hh

//---------------------------------------------------------------------------//
//                 end of PreconditionerBuilder.t.hh
//---------------------------------------------------------------------------//
