//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/Eigenvalue_Solver.t.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Mon Mar 10 14:20:33 2014
 * \brief  Eigenvalue_Solver template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_spn_Eigenvalue_Solver_t_hh
#define SPn_spn_Eigenvalue_Solver_t_hh

#include <cmath>
#include <string>

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "comm/P_Stream.hh"
#include "comm/global.hh"
#include "utils/String_Functions.hh"
#include "solvers/EigenvalueSolverBuilder.hh"
#include "solvers/InverseOperator.hh"
#include "solvers/PreconditionerBuilder.hh"
#include "Linear_System_FV.hh"
#include "Energy_Multigrid.hh"
#include "Eigenvalue_Solver.hh"
#include "MatrixTraits.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class T>
Eigenvalue_Solver<T>::Eigenvalue_Solver(RCP_ParameterList db)
    : Base(db)
    , d_keff(2.0)
{
    REQUIRE(!b_db.is_null());
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Setup the solver.
 *
 * Calls to this function builds the linear SPN system.
 */
template <class T>
void Eigenvalue_Solver<T>::setup(RCP_Dimensions  dim,
                                 RCP_Mat_DB      mat,
                                 RCP_Mesh        mesh,
                                 RCP_Indexer     indexer,
                                 RCP_Global_Data data,
                                 bool            adjoint)
{
    REQUIRE(!b_db.is_null());
    REQUIRE(!dim.is_null());
    REQUIRE(!mat.is_null());
    REQUIRE(!mesh.is_null());
    REQUIRE(!indexer.is_null());
    REQUIRE(!data.is_null());

    d_mat = mat;
    set_default_parameters();
    REQUIRE(b_db->isSublist("eigenvalue_db"));

    // build the linear system (we only provide finite volume for now)
    std::string &eqn_type =
        b_db->template get<std::string>("eqn_type", std::string("fv"));

    if (profugus::to_lower(eqn_type) == "fv")
    {
        b_system = Teuchos::rcp(
            new Linear_System_FV<T>(
                b_db, dim, mat, mesh, indexer, data));
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
        d_u = VectorTraits<T>::build_vector(Base::b_system->get_Map());

        // initialize it to 1.0
        VectorTraits<T>::put_scalar(d_u,1.0);
    }
    else
    {
        // If we already have a vector, we need to temporarily store it
        //  get the new vector, then copy the data from the old vector
        //  into the new one.  We can't just use the old vector because
        //  its internal map was destroyed when the new matrices were built,
        //  which leads to problems when something downstream of here needs
        //  access to the map.
        RCP_MV tmp_vec = d_u;

        // make the eigenvector
        d_u = VectorTraits<T>::build_vector(Base::b_system->get_Map());

        // Assign data from old vector to new one
        Teuchos::ArrayRCP<const double> tmp_data =
            VectorTraits<T>::get_data(tmp_vec);
        Teuchos::ArrayRCP<double> u_data =
            VectorTraits<T>::get_data_nonconst(d_u);
        u_data.deepCopy(tmp_data());
    }

    CHECK(!d_u.is_null());

    // get the eigenvalue solver settings
    RCP_ParameterList edb = Teuchos::sublist(b_db, "eigenvalue_db");

    // set adjoint
    b_system->set_adjoint(adjoint);

    // Build a preconditioenr
    RCP_OP prec = build_preconditioner(dim, mat, mesh, indexer, data);

    // Build the eigensolver
    d_eigensolver = EigenvalueSolverBuilder<T>::build_solver(
        edb, b_system->get_Operator(), b_system->get_fission_matrix(), prec);

    ENSURE(!d_eigensolver.is_null());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Setup the solver from an existing linear system.
 *
 */
template <class T>
void Eigenvalue_Solver<T>::setup(RCP_Mat_DB        mat,
                                 RCP_Mesh          mesh,
                                 RCP_Indexer       indexer,
                                 RCP_Global_Data   data,
                                 RCP_Linear_System system,
                                 bool              adjoint)
{
    REQUIRE(!b_db.is_null());
    REQUIRE(!mat.is_null());
    REQUIRE(!mesh.is_null());
    REQUIRE(!indexer.is_null());
    REQUIRE(!data.is_null());
    REQUIRE(!system.is_null());
    REQUIRE(!system->get_Operator().is_null());
    REQUIRE(!system->get_fission_matrix().is_null());

    // assign the material
    d_mat = mat;

    // setup default parameters
    set_default_parameters();
    CHECK(b_db->isSublist("eigenvalue_db"));

    // assign the linear system
    Base::b_system = system;

    // If we haven't built an eigenvector yet, build one and initialize it
    if( d_u.is_null() )
    {
        // make the eigenvector
        d_u = VectorTraits<T>::build_vector(Base::b_system->get_Map());

        // initialize it to 1.0
        VectorTraits<T>::put_scalar(d_u,1.0);
    }

    CHECK(!d_u.is_null());

    // get the eigenvalue solver settings
    RCP_ParameterList edb = Teuchos::sublist(b_db, "eigenvalue_db");

    // set adjoint
    b_system->set_adjoint(adjoint);

    // Build a preconditioenr
    RCP_OP prec = build_preconditioner(
        b_system->get_dims(), mat, mesh, indexer, data);

    // Build the eigensolver
    d_eigensolver = EigenvalueSolverBuilder<T>::build_solver(
        edb, b_system->get_Operator(), b_system->get_fission_matrix(), prec);

    ENSURE(!d_eigensolver.is_null());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the SPN eigenvalue equations.
 */
template <class T>
void Eigenvalue_Solver<T>::solve(Teuchos::RCP<const External_Source> q)
{
    REQUIRE( q == Teuchos::null );
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
template <class T>
void Eigenvalue_Solver<T>::write_state(State &state)
{
    REQUIRE(state.mesh().num_cells() *
             b_system->get_dims()->num_equations() * state.num_groups()
             <= VectorTraits<T>::local_length(d_u) );

    Base::write_u_into_state(d_u, state);

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
        State::View_Field moments = state.flux(g, g);

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
        State::View_Field moments = state.flux(g, g);

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
 * \brief Set default db entries for eigenvalue solvers
 */
template <>
void Eigenvalue_Solver<EpetraTypes>::set_default_parameters()
{
    // Get eigenvalue db
    RCP_ParameterList eig_db = Teuchos::sublist(
        b_db, std::string("eigenvalue_db"));
    CHECK(!eig_db.is_null());

    // The default energy multigrid preconditioner fails on one group.
    // Switch to something else.
    std::string prec_type = profugus::to_lower(
        eig_db->get("Preconditioner",std::string("Multigrid")));
    if ( d_mat->xs().num_groups() == 1 &&
         (prec_type=="multigrid" || prec_type=="multilevel") )
    {
        eig_db->set("Preconditioner", std::string("Ifpack"));
    }

    // Look for user specified tolerance in a few places, set a default if we
    // can't find one
    double tol = b_db->get("tolerance", 1.0e-6);
    tol = eig_db->get("tolerance", tol);
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

template <>
void Eigenvalue_Solver<TpetraTypes>::set_default_parameters()
{
    // Get eigenvalue db
    RCP_ParameterList eig_db = Teuchos::sublist(
        b_db, std::string("eigenvalue_db"));
    CHECK(!eig_db.is_null());

    // The default energy multigrid preconditioner fails on one group.
    // Switch to something else.
    std::string prec_type = profugus::to_lower(
        eig_db->get("Preconditioner",std::string("Ifpack2")));
    if ( d_mat->xs().num_groups() == 1 &&
         (prec_type=="multigrid" || prec_type=="multilevel") )
    {
        eig_db->set("Preconditioner", std::string("Ifpack2"));
    }

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
    // If no parameters are provided, the default solver is a MueLu
    //  algebraic multigrid preconditioner.
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
  <Parameter name='Preconditioner' type='string' value='Ifpack2'/>          \n\
  <ParameterList name='Multigrid Preconditioner'>                           \n\
   <ParameterList name='Smoother'>                                          \n\
    <Parameter name='Preconditioner' type='string' value='Ifpack2'/>        \n\
    <Parameter name='Ifpack2_Type' type='string' value='ILUT'/>             \n\
    <Parameter name='verbosity' type='string' value='None'/>                \n\
    <Parameter name='max_itr' type='int' value='5'/>                        \n\
    <Parameter name='solver_type' type='string' value='stratimikos'/>       \n\
    <ParameterList name='Stratimikos'>                                      \n\
     <Parameter name='Linear Solver Type' type='string' value='Belos'/>     \n\
     <Parameter name='Preconditioner Type' type='string' value='None'/>     \n\
    </ParameterList>                                                        \n\
   <ParameterList name='Ifpack2 Params'>                                    \n\
    <Parameter name='fact: drop tolerance' type='double' value='1e-1'/>     \n\
   </ParameterList>                                                         \n\
   </ParameterList>                                                         \n\
  </ParameterList>                                                          \n\
  <ParameterList name='Anasazi'>                                            \n\
    <Parameter name='Maximum Subspace Dimension' type='int' value='25'/>    \n\
    <Parameter name='Restart Dimension' type='int' value='5'/>              \n\
    <Parameter name='Output Level' type='string' value='medium'/>           \n\
  </ParameterList>                                                          \n\
  <Parameter name='Ifpack2_Type' type='string' value='ILUT'/>               \n\
  <ParameterList name='Ifpack2 Params'>                                     \n\
   <Parameter name='fact: drop tolerance' type='double' value='1e-2'/>      \n\
   <Parameter name='fact: ilut level-of-fill' type='double' value='1.2'/>   \n\
  </ParameterList>                                                          \n\
  <ParameterList name='MueLu Params'>                                       \n\
   <Parameter name='multigrid algorithm' type='string' value='unsmoothed'/> \n\
   <Parameter name='max levels' type='int' value='5'/>                      \n\
   <Parameter name='print initial parameters' type='bool' value='false'/>   \n\
   <Parameter name='verbosity' type='string' value='low'/>                  \n\
   <Parameter name='smoother: type' type='string' value='ILUT'/>            \n\
   <ParameterList name='smoother: params'>                                  \n\
    <Parameter name='relaxation: type' type='string' value='Gauss-Seidel'/> \n\
    <Parameter name='relaxation: sweeps' type='int' value='5'/>             \n\
   </ParameterList>                                                         \n\
   <Parameter name='coarse: type' type='string' value='ILUT'/>              \n\
   <ParameterList name='coarse: params'>                                    \n\
    <Parameter name='relaxation: type' type='string' value='Gauss-Seidel'/> \n\
    <Parameter name='relaxation: sweeps' type='int' value='5'/>             \n\
   </ParameterList>                                                         \n\
  </ParameterList>                                                          \n\
  <ParameterList name='operator_db'>                                        \n\
   <Parameter name='verbosity' type='string' value='Low'/>                  \n\
   <Parameter name='solver_type' type='string' value='stratimikos'/>        \n\
   <ParameterList name='Stratimikos'>                                       \n\
    <Parameter name='Linear Solver Type' type='string' value='Belos'/>      \n\
    <Parameter name='Preconditioner Type' type='string' value='None'/>      \n\
   </ParameterList>                                                         \n\
  </ParameterList>                                                          \n\
 </ParameterList>                                                           \n"
        );

    // Convert string to a Teuchos PL
    RCP_ParameterList default_pl =
        Teuchos::getParametersFromXmlString(plrefstr);

    // Insert defaults into pl, leaving existing values in tact
    eig_db->setParametersNotAlreadySet(*default_pl);

    // If Ifpack is selected, change it to Ifpack2
    if( eig_db->get<std::string>("Preconditioner")=="Ifpack" )
    {
        eig_db->set("Preconditioner",std::string("Ifpack2"));
    }

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
template <class T>
Teuchos::RCP<typename T::OP>
Eigenvalue_Solver<T>::build_preconditioner(RCP_Dimensions  dim,
                                           RCP_Mat_DB      mat,
                                           RCP_Mesh        mesh,
                                           RCP_Indexer     indexer,
                                           RCP_Global_Data data)
{
    REQUIRE(b_db->isSublist("eigenvalue_db"));

    // preconditioner operator
    RCP_OP prec;

    // get the eigenvalue database
    RCP_ParameterList edb = Teuchos::sublist(b_db, "eigenvalue_db");

    // get the preconditioner type
    std::string prec_type = profugus::to_lower(
        edb->template get<std::string>("Preconditioner",
                                       std::string("Multigrid")));

    // build the appropriate preconditioners
    if (prec_type=="multigrid" || prec_type=="multilevel")
    {
        RCP_ParameterList prec_db = Teuchos::sublist(
            edb, "Multigrid Preconditioner");

        prec = Teuchos::rcp(
            new Energy_Multigrid<T>(b_db, prec_db, dim, mat, mesh,
                                    indexer, data, b_system));
        CHECK(prec != Teuchos::null);
    }
    else if (prec_type == "stratimikos" || prec_type == "solver")
    {
        // Create default db
        RCP_ParameterList op_db = Teuchos::sublist(edb, "operator_db");

        // Build Stratimikos Operator
        Teuchos::RCP<InverseOperator<T> > solver_op = Teuchos::rcp(
                new InverseOperator<T>(op_db));
        solver_op->set_operator(b_system->get_Operator());

        prec = solver_op;

        CHECK(prec != Teuchos::null);
    }
    else
    {
        prec = PreconditionerBuilder<T>::build_preconditioner(
            b_system->get_Operator(),edb);
    }

    return prec;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write matrices to file
 */
template <class T>
void Eigenvalue_Solver<T>::write_problem_to_file() const
{
    // Write A (LHS operator)
    Teuchos::RCP<const typename T::MATRIX> matrix =
        Teuchos::rcp_dynamic_cast<const typename T::MATRIX>(
            b_system->get_Operator());
    MatrixTraits<T>::write_matrix_file(matrix,"A.mtx");

    // Write B (RHS operator)
    matrix = Teuchos::rcp_dynamic_cast<const typename T::MATRIX>(
        b_system->get_fission_matrix());
    MatrixTraits<T>::write_matrix_file(matrix,"B.mtx");
}

} // end namespace profugus

#endif // SPn_spn_Eigenvalue_Solver_t_hh

//---------------------------------------------------------------------------//
// end of Eigenvalue_Solver.t.hh
//---------------------------------------------------------------------------//
