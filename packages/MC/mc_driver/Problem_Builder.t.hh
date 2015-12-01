//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_driver/Problem_Builder.t.hh
 * \author Thomas M. Evans
 * \date   Wed Mar 12 22:25:22 2014
 * \brief  Problem_Builder member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_driver_Problem_Builder_t_hh
#define mc_driver_Problem_Builder_t_hh

#include <utility>
#include <algorithm>
#include <sstream>

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "utils/String_Functions.hh"
#include "xs/XS_Builder.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "geometry/Mesh_Geometry.hh"
#include "mc/Box_Shape.hh"
#include "mc/VR_Analog.hh"
#include "mc/VR_Roulette.hh"
#include "mc/Fission_Matrix_Tally.hh"
#include "mc/Source_Diagnostic_Tally.hh"
#include "Problem_Builder.hh"
#include "Geometry_Builder.hh"

namespace mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Problem_Builder<Geometry>::Problem_Builder()
{
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Setup the problem.
 */
template <class Geometry>
void Problem_Builder<Geometry>::setup(RCP_ParameterList master)
{
    // validate the parameter list
    INSIST(master->isSublist("MATERIAL"),
            "MATERIAL block not defined in input.");
    INSIST(master->isSublist("PROBLEM"),
            "PROBLEM block not defined in input.");

    // store the individual parameter lists
    d_matdb    = Teuchos::sublist(master, "MATERIAL");
    d_db       = Teuchos::sublist(master, "PROBLEM");

    CHECK(!d_db.is_null());

    // default the boundary conditions to reflecting
    if (!d_db->isParameter("boundary"))
    {
        d_db->set("boundary", "reflect");
    }

    // add the default boundary treatment for reflecting
    if (d_db->get<std::string>("boundary") == "reflect")
    {
        if (!d_db->isSublist("boundary_db"))
        {
            // set to all reflecting b.c. if undefined
            Teuchos::ParameterList boundary("boundary");
            OneDArray_int reflect(6, 1);
            boundary.set("reflect", reflect);

            d_db->set("boundary_db", boundary);
        }
    }

    // build the geometry
    Geometry_Builder<Geometry> geom_builder;
    d_geometry = geom_builder.build(master);

    // build build physics
    build_physics();

    // build the variance reduction
    build_var_reduction();

    // build the external source (there won't be one for k-eigenvalue
    // problems)
    if (master->isSublist("SOURCE"))
    {
        build_source(master->sublist("SOURCE"));
    }

    // build the SPN problem for fission matrix acceleration
    build_spn_problem();

    // build the tallier
    build_tallies();
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Build the Monte Carlo MG physics.
 *
 * For now, we build all the cross sections in the problem on every domain.
 * For the mini-app, this is not expected to be an overburdening cost.
 */
template <class Geometry>
void Problem_Builder<Geometry>::build_physics()
{
    typedef profugus::XS_Builder::Matid_Map Matid_Map;
    typedef profugus::XS_Builder::RCP_XS    RCP_XS;

    REQUIRE(d_matdb->isParameter("mat list"));
    VALIDATE(d_matdb->isParameter("xs library"),
              "Inline cross sections not implemented yet.");

    // get the material list off of the database
    const auto &mat_list = d_matdb->get<OneDArray_str>("mat list");

    // convert the matlist to a mat-id map
    Matid_Map matids;
    for (int id = 0, N = mat_list.size(); id < N; ++id)
    {
        matids.insert(Matid_Map::value_type(id, mat_list[id]));
    }
    matids.complete();
    CHECK(matids.size() == mat_list.size());

    // make a cross section builder
    profugus::XS_Builder builder;

    // broadcast the raw cross section data
    builder.open_and_broadcast(d_matdb->get<std::string>("xs library"));

    // get the number of groups and moments in the cross section data
    int Ng_data = builder.num_groups();
    int N_data  = builder.pn_order();

    // get the number of groups required
    int g_first = d_db->get("g_first", 0);
    int g_last  = d_db->get("g_last", Ng_data - 1);
    VALIDATE(1 + (g_last - g_first) <= Ng_data, "Energy group range exceeds "
              << "number of groups in data, 1 + g_last - g_first = "
              << 1 + (g_last - g_first) << " > " << Ng_data);

    // build the cross sections (always build P0 for Monte Carlo)
    builder.build(matids, 0, g_first, g_last);
    RCP_XS xs = builder.get_xs();
    CHECK(xs->num_mat() == matids.size());
    CHECK(xs->num_groups() == 1 + (g_last - g_first));

    // make the physics
    d_physics = std::make_shared<Physics_t>(d_db, xs);

    // set the geometry in the physics
    d_physics->set_geometry(d_geometry);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the variance reduction.
 */
template <class Geometry>
void Problem_Builder<Geometry>::build_var_reduction()
{
    using profugus::to_lower;

    REQUIRE(!d_db.is_null());

    // the default is to do roulette
    const auto &var = to_lower(
        d_db->get<std::string>("variance reduction", std::string("roulette")));

    // build the appropriate variance reduction
    if (to_lower(var) == "roulette")
    {
        d_var_reduction =
            std::make_shared<profugus::VR_Roulette<Geom_t> >(d_db);
    }
    else if (to_lower(var) == "analog")
    {
        d_var_reduction =
            std::make_shared<profugus::VR_Analog<Geom_t> >();
    }
    else
    {
        std::ostringstream m;
        m << "Variance reduction type " << var << "unknown.";
        throw profugus::assertion(m.str());
    }

    ENSURE(d_var_reduction);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the external source.
 *
 *\param source_db
 */
template <class Geometry>
void Problem_Builder<Geometry>::build_source(const ParameterList &source_db)
{
    REQUIRE(source_db.isParameter("box"));
    REQUIRE(source_db.isParameter("spectrum"));
    REQUIRE(d_physics);

    // get the source box coordinates
    const auto &box = source_db.get<OneDArray_dbl>("box");
    CHECK(box.size() == 6);

    // make the box shape
    d_shape = std::make_shared<profugus::Box_Shape>(
        box[0], box[1], box[2], box[3], box[4], box[5]);

    // get the source spectrum and add it to the main database
    const auto &shape = source_db.get<OneDArray_dbl>("spectrum");
    CHECK(shape.size() == d_physics->num_groups());

    // add the shape to the main database because the MC source gets the
    // spectral shape from the main db
    d_db->set("spectral_shape", shape);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the tallies.
 */
template <class Geometry>
void Problem_Builder<Geometry>::build_tallies()
{
    using profugus::Mesh_Geometry;

    typedef profugus::Fission_Matrix_Tally<Geom_t>    FM_Tally_t;
    typedef profugus::Source_Diagnostic_Tally<Geom_t> SD_Tally_t;

    // make the tallier
    d_tallier = std::make_shared<Tallier_t>();
    d_tallier->set(d_geometry, d_physics);

    // check for fission matrix tallies
    if (d_db->isSublist("fission_matrix_db"))
    {
        // get the database
        const ParameterList &fdb = d_db->sublist("fission_matrix_db");

        // if this is a tally then add it
        if (fdb.isSublist("tally"))
        {
            INSIST(!d_shape,
                   "Cannot run fission matrix on fixed-source problems.");

            // get the database
            const ParameterList &tdb = fdb.sublist("tally");

            // validate that there is a fission matrix mesh defined
            VALIDATE(tdb.isParameter("x_bounds"),
                     "Failed to define x-boundaries for fission matrix mesh");
            VALIDATE(tdb.isParameter("y_bounds"),
                     "Failed to define y-boundaries for fission matrix mesh");
            VALIDATE(tdb.isParameter("z_bounds"),
                     "Failed to define z-boundaries for fission matrix mesh");

            // get the fission matrix mesh boundaries
            auto xb = tdb.get<OneDArray_dbl>("x_bounds").toVector();
            auto yb = tdb.get<OneDArray_dbl>("y_bounds").toVector();
            auto zb = tdb.get<OneDArray_dbl>("z_bounds").toVector();
            typename FM_Tally_t::SP_Mesh_Geometry geo(
                std::make_shared<Mesh_Geometry>(xb, yb, zb));

            // build the fission matrix
            auto fm_tally(std::make_shared<FM_Tally_t>(
                              d_db, d_physics, geo));
            CHECK(fm_tally);

            // add this to the tallier
            d_tallier->add_compound_tally(fm_tally);
        }
    }

    // check for source tallies
    if (d_db->isSublist("source_diagnostic_db"))
    {
        INSIST(!d_shape,
               "Cannot run source diagnostic on fixed-source problems.");

        // get the database
        const ParameterList &sdb = d_db->sublist("source_diagnostic_db");

        // validate that there is a source diagnostic mesh defined
        VALIDATE(sdb.isParameter("x_bounds"),
                 "Failed to define x-boundaries for source diagnostic mesh");
        VALIDATE(sdb.isParameter("y_bounds"),
                 "Failed to define y-boundaries for source diagnostic mesh");
        VALIDATE(sdb.isParameter("z_bounds"),
                 "Failed to define z-boundaries for source diagnostic mesh");

        // get the fission matrix mesh boundaries
        auto xb = sdb.get<OneDArray_dbl>("x_bounds").toVector();
        auto yb = sdb.get<OneDArray_dbl>("y_bounds").toVector();
        auto zb = sdb.get<OneDArray_dbl>("z_bounds").toVector();
        typename SD_Tally_t::SP_Mesh_Geometry geo(
            std::make_shared<Mesh_Geometry>(xb, yb, zb));

        // the default is to tally during inactive cycles

        // build the source tally
        auto src_tally(std::make_shared<SD_Tally_t>(
                           d_db, d_physics, geo, true));
        CHECK(src_tally);

        // add this to the tallier
        d_tallier->add_source_tally(src_tally);
    }

    ENSURE(d_tallier);
}

//---------------------------------------------------------------------------//
/*!
 * \build the SPN problem
 */
template <class Geometry>
void Problem_Builder<Geometry>::build_spn_problem()
{
    typedef profugus::EpetraTypes                          ET;
    typedef profugus::TpetraTypes                          TT;
    typedef profugus::Fission_Matrix_Acceleration_Impl<Geom_t,ET> FM_ET;
    typedef profugus::Fission_Matrix_Acceleration_Impl<Geom_t,TT> FM_TT;

    // check for fission matrix acceleration
    if (d_db->isSublist("fission_matrix_db"))
    {
        // get the database
        const ParameterList &fdb = d_db->sublist("fission_matrix_db");

        // if this is an acceleration, then build the SPN problem builder
        if (fdb.isSublist("acceleration"))
        {
            INSIST(!d_shape,
                   "Cannot run fission matrix on fixed-source problems.");

            // get the database
            const ParameterList &adb = fdb.sublist("acceleration");

            // validate parameters
            VALIDATE(adb.isParameter("spn_problem"),
                     "Failed to define the equivalent SPN problem for "
                     << "fission matrix acceleration.");

            // get the xml-file for the SPN problem
            auto spn_input = adb.get<std::string>("spn_problem");

            // make the SPN problem builder
            SPN_Builder spn_builder;

            // setup the SPN problem
            spn_builder.setup(spn_input);

            // get the SPN problem database
            auto spn_db = spn_builder.problem_db();

            // get the linear algebra type
            std::string type = spn_db->get("trilinos_implementation",
                                           std::string("epetra"));

            // build the fission matrix acceleration
            if (type == "epetra")
            {
                d_fm_acceleration = std::make_shared<FM_ET>();
            }
            else if (type == "tpetra")
            {
                d_fm_acceleration = std::make_shared<FM_TT>();
            }

            // build the problem
            d_fm_acceleration->build_problem(spn_builder);

            // get the global mesh
            const auto &mesh = d_fm_acceleration->global_mesh();

            // make the bounds
            OneDArray_dbl x(mesh.edges(0));
            OneDArray_dbl y(mesh.edges(1));
            OneDArray_dbl z(mesh.edges(2));

            // make a source tally that defaults to the acceleration mesh if a
            // source tally hasn't been defined
            ParameterList &sdb = d_db->sublist("source_diagnostic_db");

            // make a default parameterlist
            ParameterList def;
            def.set("x_bounds", x);
            def.set("y_bounds", y);
            def.set("z_bounds", z);

            // add these to the source diagnosic tally
            sdb.setParametersNotAlreadySet(def);
        }
    }
}

} // end namespace mc

#endif // mc_driver_Problem_Builder_t_hh

//---------------------------------------------------------------------------//
//                 end of Problem_Builder.t.hh
//---------------------------------------------------------------------------//
