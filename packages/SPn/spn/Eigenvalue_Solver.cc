//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Eigenvalue_Solver.cc
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Mon Mar 10 14:20:33 2014
 * \brief  Eigenvalue_Solver member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <string>

#include "Ifpack.h"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Epetra_InvOperator.h"

#include "harness/DBC.hh"
#include "comm/P_Stream.hh"
#include "comm/global.hh"
#include "utils/String_Functions.hh"
#include "solvers/EigenvalueSolverBuilder.hh"
#include "solvers/InverseOperator.hh"
#include "Linear_System_FV.hh"
#include "Energy_Multigrid.hh"
#include "Eigenvalue_Solver.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Eigenvalue_Solver::Eigenvalue_Solver(RCP_ParameterList db)
    : Base(db)
    , d_keff(2.0)
{
    REQUIRE(!b_db.is_null());

    set_default_parameters();
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Setup the solver.
 *
 * Calls to this function builds the linear SPN system.
 */
void Eigenvalue_Solver::setup(RCP_Dimensions  dim,
                              RCP_Mat_DB      mat,
                              RCP_Mesh        mesh,
                              RCP_Indexer     indexer,
                              RCP_Global_Data data)
{
    REQUIRE(!b_db.is_null());
    REQUIRE(!dim.is_null());
    REQUIRE(!mat.is_null());
    REQUIRE(!mesh.is_null());
    REQUIRE(!indexer.is_null());
    REQUIRE(!data.is_null());
    REQUIRE(b_db->isSublist("eigenvalue_db"));

    // build the linear system (we only provide finite volume for now)
    std::string &eqn_type = b_db->get("eqn_type", std::string("fv"));

    if (profugus::to_lower(eqn_type) == "fv")
    {
        b_system = Teuchos::rcp(
            new Linear_System_FV(b_db, dim, mat, mesh, indexer, data));
    }
    else
    {
        std::string msg = "Undefined equation type: " + eqn_type;
        throw profugus::assertion(msg);
    }
    CHECK(!b_system.is_null());

    // build the SPN matrix and build the fission matrix
    b_system->build_Matrix();
    b_system->build_fission_matrix();

    // If we haven't built an eigenvector yet, build one and initialize it
    // If we already have one, copy its data into new vector
    if( d_u.is_null() )
    {
        // make the eigenvector
        d_u = Teuchos::rcp(new Vector_t(*Base::b_system->get_Map()));

        // initialize it to 1.0
        d_u->PutScalar(1.0);
    }
    else
    {
        // If we already have a vector, we need to temporarily store it
        //  get the new vector, then copy the data from the old vector
        //  into the new one.  We can't just use the old vector because
        //  its internal map was destroyed when the new matrices were built,
        //  which leads to problems when something downstream of here needs
        //  access to the map.
        RCP_Vector tmp_vec = d_u;

        // make the eigenvector
        d_u = Teuchos::rcp(new Vector_t(*Base::b_system->get_Map()));

        // Assign data from old vector to new one
        d_u->Update(1.0, *tmp_vec, 0.0);
    }

    // the map is a block map where the block size is the number of groups
    //  OR a point map
    CHECK(b_system->get_Map()->NumMyElements() * mat->xs().num_groups()
           == d_u->MyLength() ||
           b_system->get_Map()->NumMyElements() == d_u->MyLength() );
    CHECK(!d_u.is_null());

    // get the eigenvalue solver settings
    RCP_ParameterList edb = Teuchos::sublist(b_db, "eigenvalue_db");

    // The default energy multigrid preconditioner for the Davidson solver
    // fails on one group.  We could switch to a different preconditioner but
    // we would have to set up parameters for it.  Instead, we'll just switch
    // to Arnoldi for that case.
    if ( mat->xs().num_groups() == 1 )
    {
        edb->set("Preconditioner", std::string("Ifpack"));
    }

    // Build a preconditioenr
    RCP_Epetra_Op prec = build_preconditioner(dim, mat, mesh, indexer, data);

    // Build the eigensolver
    d_eigensolver = EigenvalueSolverBuilder<EPETRA>::build_solver(
        edb, b_system->get_Operator(), b_system->get_fission_matrix(), prec);

    ENSURE(!d_eigensolver.is_null());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the SPN eigenvalue equations.
 */
void Eigenvalue_Solver::solve()
{
    REQUIRE(!b_system.is_null());
    REQUIRE(!d_u.is_null());
    REQUIRE(!d_eigensolver.is_null());

    // solve the problem
    d_eigensolver->solve(d_keff, d_u);

    profugus::pout << profugus::setprecision(10) << profugus::fixed;
    profugus::pout << "k-eff = " << d_keff << profugus::endl;
    profugus::pout << profugus::setprecision(6);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write scalar flux (eigenvector) into the state.
 */
void Eigenvalue_Solver::write_state(State_t &state)
{
    REQUIRE(state.mesh().num_cells() *
             b_system->get_dims()->num_equations() * state.num_groups()
             <= d_u->MyLength());

    Base::write_u_into_state(*d_u, state);

    // normalize the eigenvector (over all groups/space)
    double norm[2] = {0.0, 0.0};

    // flux in state at each moment/unknown/cell/energy
    double flux = 0.0;

    // number of groups
    int Ng = state.num_groups();

    // number of cells in this domain
    int Nc = state.mesh().num_cells();

    // loop over groups
    for (int g = 0; g < Ng; ++g)
    {
        // get moments for this group
        View_Field moments = state.flux(g, g);

        // loop over cells on this domain
        for (int cell = 0; cell < Nc; ++cell)
        {
            // get the scalar flux at this point
            flux = moments[cell];

            // add to the normalization factor
            norm[0] += flux * flux;

            // sum the vector to determine the sign of the eigenvector
            norm[1] += flux;
        }
    }

    // do a global reduction on the normalization
    profugus::global_sum(norm, 2);
    CHECK(norm[0] > 0.0);

    // apply the normalization
    double norm_f = (1.0 / std::sqrt(norm[0])) * (std::fabs(norm[1]) / norm[1]);

    // loop over groups
    for (int g = 0; g < Ng; ++g)
    {
        // get moments for this group
        View_Field moments = state.flux(g, g);

        // loop over cells on this domain
        for (int cell = 0; cell < Nc; ++cell)
        {
            // get the scalar flux at this point
            moments[cell] *= norm_f;
        }
    }
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set default db entries for Davidson solver
 */
void Eigenvalue_Solver::set_default_parameters()
{
    // Get eigenvalue db
    RCP_ParameterList eig_db = Teuchos::sublist(
        b_db, std::string("eigenvalue_db"));
    CHECK(!eig_db.is_null());

    // Look for user specified tolerance in a few places, set a default if we
    // can't find one
    double tol = b_db->get("tolerance", 1.0e-6);
    tol        = eig_db->get("tolerance", tol);
    eig_db->get("Convergence Tolerance", tol);

    // look for user specified iteration limit in a few places, set a default
    // if we can't find one
    int max_itr = b_db->get("max_itr", 500);
    max_itr     = eig_db->get("max_itr", max_itr);
    eig_db->get("Maximum Restarts", max_itr);

    // Forgive the formatting, column limits are a bit of an issue
    // If no parameters are provided, the default solver is a
    //  Davidson solver with a multigrid-in-energy preconditioner which
    //  internally uses an ilu-preconditioned BiCGStab as a smoother
    // The goal here is to support "good" default parameters for a variety
    //  of different solvers and preconditioners so that the user can toggle
    //  between a small number of options without having to set many
    //  parameters.  Advanced users can still specify as much detail
    //  as they want.
    const std::string plrefstr(
        "<ParameterList name='eigenvalue_db'>                               \n\
  <Parameter name='eigensolver' type='string' value='Davidson'/>            \n\
  <Parameter name='Output Level' type='string' value='low'/>                \n\
  <Parameter name='verbosity' type='string' value='Low'/>                   \n\
  <Parameter name='Preconditioner' type='string' value='Multigrid'/>        \n\
  <ParameterList name='Multigrid Preconditioner'>                           \n\
   <ParameterList name='Smoother'>                                          \n\
    <Parameter name='Preconditioner' type='string' value='Ifpack'/>         \n\
    <Parameter name='verbosity' type='string' value='None'/>                \n\
    <Parameter name='max_itr' type='int' value='3'/>                        \n\
    <Parameter name='aztec_solver' type='string' value='BiCGStab'/>         \n\
    <Parameter name='aztec_prec' type='string' value='ilu'/>                \n\
   </ParameterList>                                                         \n\
  </ParameterList>                                                          \n\
  <ParameterList name='Anasazi'>                                            \n\
    <Parameter name='Maximum Subspace Dimension' type='int' value='25'/>    \n\
    <Parameter name='Restart Dimension' type='int' value='5'/>              \n\
  </ParameterList>                                                          \n\
  <Parameter name='Ifpack Type' type='string' value='ILUT'/>                \n\
  <Parameter name='Ifpack Overlap' type='int' value='0'/>                   \n\
  <ParameterList name='Ifpack Params'>                                      \n\
    <Parameter name='fact: drop tolerance' type='double' value='1e-2'/>     \n\
    <Parameter name='fact: ilut level-of-fill' type='double' value='1.2'/>  \n\
  </ParameterList>                                                          \n\
  <Parameter name='ML Default Type' type='string' value='DD'/>              \n\
  <ParameterList name='ML Params'>                                          \n\
    <Parameter name='smoother: type' type='string' value='ILU'/>            \n\
    <Parameter name='smoother: sweeps' type='int' value='3'/>               \n\
    <Parameter name='smoother: ifpack overlap' type='int' value='0'/>       \n\
    <Parameter name='max levels' type='int' value='4'/>                     \n\
  </ParameterList>                                                          \n\
  <ParameterList name='operator_db'>                                        \n\
   <Parameter name='verbosity' type='string' value='Low'/>                  \n\
   <Parameter name='aztec_solver' type='string' value='bicgstab'/>          \n\
   <Parameter name='aztec_prec' type='string' value='ilu'/>                 \n\
  </ParameterList>                                                          \n\
 </ParameterList>                                                           \n"
        );

    // Convert string to a Teuchos PL
    RCP_ParameterList default_pl =
        Teuchos::getParametersFromXmlString(plrefstr);

    // Insert defaults into pl, leaving existing values in tact
    eig_db->setParametersNotAlreadySet(*default_pl);

    // propagate stopping criteria for Arnoldi
    eig_db->sublist("Anasazi").get("Convergence Tolerance", tol);
    eig_db->sublist("Anasazi").get("Maximum Restarts", max_itr);

    // propagate stopping criteria for operators
    eig_db->sublist("operator_db").get("tolerance", 0.1 * tol);
    eig_db->sublist("operator_db").get("max_itr", max_itr);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build preconditioner
 */
Eigenvalue_Solver::RCP_Epetra_Op
Eigenvalue_Solver::build_preconditioner(RCP_Dimensions  dim,
                                        RCP_Mat_DB      mat,
                                        RCP_Mesh        mesh,
                                        RCP_Indexer     indexer,
                                        RCP_Global_Data data)
{
    REQUIRE(b_db->isSublist("eigenvalue_db"));

    // preconditioner operator
    RCP_Epetra_Op prec;

    // get the eigenvalue database
    RCP_ParameterList edb = Teuchos::sublist(b_db, "eigenvalue_db");

    // get the preconditioner type
    std::string prec_type = profugus::to_lower(
        edb->get("Preconditioner", std::string("Multigrid")));

    // build the appropriate preconditioners
    if (prec_type=="multigrid" || prec_type=="multilevel")
    {
        RCP_ParameterList prec_db = Teuchos::sublist(
            edb, "Multigrid Preconditioner");

        prec = Teuchos::rcp(
            new Energy_Multigrid(b_db, prec_db, dim, mat, mesh,
                                 indexer, data, b_system));
        CHECK(prec != Teuchos::null);
    }
    else if (prec_type == "ifpack")
    {
        // Get the matrix from the linear system (may not be full system
        // matrix)
        Teuchos::RCP<Epetra_RowMatrix> rcp_rowmat = b_system->get_Matrix();
        INSIST(!rcp_rowmat.is_null(),
               "Cannot use Ifpack preconditioner without constructing matrix.");

        // Create Ifpack preconditioner
        Ifpack ifpack_factory;
        std::string ifpack_type = edb->get<std::string>("Ifpack Type");
        int overlap             = edb->get<int>("Ifpack Overlap");
        d_ifpack_prec           = Teuchos::rcp(
            ifpack_factory.Create(ifpack_type, rcp_rowmat.get(), overlap));

        // Convert "Ifpack Params" database to Teuchos::ParameterList
        Teuchos::ParameterList &ifpack_db = edb->sublist("Ifpack Params");
        d_ifpack_prec->SetParameters(ifpack_db);

        // Process preconditioner
        d_ifpack_prec->Initialize();
        d_ifpack_prec->Compute();
        if( edb->get<std::string>("Output Level") == "high" &&
            profugus::node() == 0 )
        {
            std::cout << "Ifpack Parameter List" << std::endl;
            d_ifpack_prec->Print(std::cout);
        }

        // Ifpack preconditioners are applied with "Apply_Inverse"
        //  but we want it to be used with "Apply", wrap the
        //  preconditioner into an object that reverses the functionality
        prec = Teuchos::rcp( new Epetra_InvOperator(d_ifpack_prec.get()) );

        CHECK( prec != Teuchos::null );
    }
    else if (prec_type == "ml")
    {
#ifdef USE_ML
        // Create default db
        RCP_ParameterList ml_db = Teuchos::sublist(edb, "ML Params");

        // Choose default settings
        std::string default_type =
            edb->get("ML Default Type", std::string("DD"));

        // Set default profile, but don't override existing entries
        std::vector<int> az_options(AZ_OPTIONS_SIZE);
        std::vector<double> az_params(AZ_PARAMS_SIZE);
        bool override = false;
        ML_Epetra::SetDefaults(default_type, *ml_db, &az_options[0],
                               &az_params[0], override);

        // Get the matrix from the linear system (may not be full system
        // matrix)
        Teuchos::RCP<Epetra_RowMatrix> rcp_rowmat = b_system->get_Matrix();
        INSIST(!rcp_rowmat.is_null(),
               "Cannot use ML preconditioner without constructing matrix.");

        d_ml_prec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(
                                     *rcp_rowmat, *ml_db));

        if (edb->get<std::string>("Output Level") == "high" &&
            profugus::node() == 0)
        {
            std::cout << "ML Defaults: " << default_type << std::endl;
            std::cout << "ML Parameter List" << std::endl;
            ml_db->print(std::cout);
        }

        // ML preconditioners are applied with "Apply_Inverse"
        //  but we want it to be used with "Apply", wrap the
        //  preconditioner into an object that reverses the functionality
        prec = Teuchos::rcp(new Epetra_InvOperator(d_ml_prec.get()));
        CHECK(!prec.is_null());
#else
        throw profugus::assertion("ML must be enabled at configure.");
#endif
    }
    else if (prec_type == "stratimikos" || prec_type == "solver")
    {
        // Create default db
        RCP_ParameterList op_db = Teuchos::sublist(edb, "operator_db");

        // Build Stratimikos Operator
        Teuchos::RCP<InverseOperator<EPETRA> > solver_op = Teuchos::rcp(
                new InverseOperator<EPETRA>(op_db));
        solver_op->set_operator(b_system->get_Operator());

        prec = solver_op;

        CHECK(prec != Teuchos::null);
    }
    else if (prec_type != "none" && prec_type != "internal")
    {
        VALIDATE(false, "Unknown preconditioner option: " << prec_type);
    }

    return prec;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Eigenvalue_Solver.cc
//---------------------------------------------------------------------------//
