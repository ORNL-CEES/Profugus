//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/PreconditionerBuilder.cc
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 09:28:52 2014
 * \brief  PreconditionerBuilder member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <string>

#include <SPn/config.h>

#include "Epetra_RowMatrix.h"
#include "Ifpack.h"

// ML has to be optional for Windows compatibility
#ifdef USE_ML
#include "ml_MultiLevelPreconditioner.h"
#endif

#include "utils/String_Functions.hh"
#include "PreconditionerBuilder.hh"

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
Teuchos::RCP<Epetra_Operator>
PreconditionerBuilder::build_preconditioner(
    Teuchos::RCP<Epetra_Operator> op,
    RCP_ParameterList             db )
{
    using std::string;

    // Default to Ifpack
    string prec_type = to_lower(db->get("Preconditioner", string("ifpack")));
    Validate(prec_type == "ifpack" || prec_type=="ml" ||
             prec_type=="none",
             "Preconditioner must be 'Ifpack', 'ML', or 'None'.");

    // Get the underlying matrix, must be Epetra_RowMatrix
    Teuchos::RCP<Epetra_RowMatrix> rowmat =
        Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(op);
    Insist( rowmat != Teuchos::null,
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
        Ensure( err == 0 );
        err = ifpack_prec->Compute();
        Ensure( err == 0 );

        // Wrap raw preconditioner into a denovo Preconditioner to
        //  invert the sense of Apply and ApplyInverse
        // This provides the same functionality as an Epetra_InvOperator
        //  but uses an RCP rather than a raw pointer so we don't
        //  have to store a copy of the raw preconditioner separately.
        prec = Teuchos::RCP<Epetra_Operator>(
            new Preconditioner(ifpack_prec) );
        Ensure( prec != Teuchos::null );
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

        // Wrap raw preconditioner into a denovo Preconditioner to
        //  invert the sense of Apply and ApplyInverse
        // This provides the same functionality as an Epetra_InvOperator
        //  but uses an RCP rather than a raw pointer so we don't
        //  have to store a copy of the raw preconditioner separately.
        prec = Teuchos::RCP<Epetra_Operator>(new Preconditioner(ml_prec) );

        Ensure( prec != Teuchos::null );
#else
        Validate(false,"ML not enabled in this build.");
#endif
    }

    return prec;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of PreconditionerBuilder.cc
//---------------------------------------------------------------------------//
