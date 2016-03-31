//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Manager_Cuda.t.cuh
 * \author Steven Hamilton
 * \date   Wed Jun 18 11:21:16 2014
 * \brief  Manager_Cuda member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Manager_Cuda_t_cuh
#define cuda_mc_Manager_Cuda_t_cuh

#include <fstream>
#include <iomanip>

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "harness/DBC.hh"
#include "comm/Timing.hh"
#include "comm/global.hh"
#include "xs/XS_Builder.hh"
#include "utils/String_Functions.hh"
#include "utils/Serial_HDF5_Writer.hh"
#include "utils/Parallel_HDF5_Writer.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "cuda_geometry/Mesh_Geometry.hh"
//#include "mc/Fission_Source.hh"
#include "Uniform_Source.cuh"
#include "Tallier.cuh"
//#include "KCode_Solver.hh"
#include "Manager_Cuda.hh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// PRIVATE TEMPLATE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Manager_Cuda<Geometry>::Manager_Cuda()
    : d_node(profugus::node())
    , d_nodes(profugus::nodes())
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Setup the problem.
 * \param master Problem parameters.
 */
template <class Geometry>
void Manager_Cuda<Geometry>::setup(RCP_ParameterList master)
{
    SCOPED_TIMER("Manager_Cuda.setup");

    SCREEN_MSG("Building and initializing geometry, physics, "
               << "variance reduction, and tallies");

    // get the problem database from the problem-builder
    d_db = Teuchos::sublist(master, "PROBLEM");
    CHECK(!d_db.is_null());

    // store the problem name
    d_problem_name = d_db->get("problem_name", std::string("MC"));
    d_db->setName(d_problem_name + "-PROBLEM");

    // get the geometry and physics
    build_geometry(master);
    CHECK(d_geometry.get_host_ptr());
    CHECK(d_geometry.get_device_ptr());

    build_physics(master);
    CHECK(d_physics.get_host_ptr());
    CHECK(d_physics.get_device_ptr());

    // get the external source shape (it could be null)
    std::string prob_type;
    if( master->isSublist("SOURCE") )
        prob_type = "fixed";
    else
        prob_type = "eigenvalue";

    // set the problem type in the final db
    d_db->set("problem_type", prob_type);

    SCREEN_MSG("Building " << prob_type << " solver");

    // get the tallier
    auto tallier_host = std::make_shared<Tallier<Geom_t>>();
    tallier_host->set(d_geometry,d_physics);

    if( d_db->isSublist("cell_tally_db") )
    {
        auto tally_db = Teuchos::sublist(d_db,"cell_tally_db");

        INSIST( tally_db->isType<Array_Int>("cells"), "Must have cell list" );

        auto cells = tally_db->get<Array_Int>("cells");

        auto cell_tally_host = std::make_shared<Cell_Tally<Geom_t>>(
            d_geometry, d_physics);
        cell_tally_host->set_cells( cells.toVector() );

        auto cell_tally =
            cuda::Shared_Device_Ptr<Cell_Tally<Geom_t>>(cell_tally_host);
        tallier_host->add_cell_tally(cell_tally);
    }

    auto tallier = cuda::Shared_Device_Ptr<Tallier<Geom_t>>(tallier_host);

    // make the transporter
    SP_Transporter transporter(std::make_shared<Transporter_t>(
                                   d_db, d_geometry, d_physics, tallier));

    // build the appropriate solver (default is eigenvalue)
    if (prob_type == "eigenvalue")
    {
#if 0
        // make the fission source
        SP_Fission_Source source(
            std::make_shared<Fission_Source_t>(
                d_db, d_geometry, d_physics, d_rnd_control));

        // make the kcode solver
        auto kcode_solver =
            std::make_shared<profugus::KCode_Solver<Geom_t> >(d_db);

        // set the solver
        kcode_solver->set(transporter, source);

        // set hybrid acceleration
        kcode_solver->set(builder.get_acceleration());

        // assign the base solver
        d_keff_solver = kcode_solver;

        // assign the base solver
        d_solver = d_keff_solver;
#else
        INSIST(false,"KCODE solver for Cuda not yet implemented.");
#endif

    }
    else if (prob_type == "fixed")
    {
        // Get source db
        auto source_db = Teuchos::sublist(master,"SOURCE");
        REQUIRE( source_db->isType<Array_Dbl>("box") );
        REQUIRE( source_db->isType<Array_Dbl>("spectrum") );

        // Source bounding box
        const auto &box = source_db->get<Array_Dbl>("box");
        CHECK( box.size() == 6 );

        d_shape = cuda::shared_device_ptr<cuda_mc::Box_Shape>(
            box[0], box[1], box[2], box[3], box[4], box[5] );

        // Source spectrum
        const auto &spectrum = source_db->get<Array_Dbl>("spectrum");
        CHECK( spectrum.size() == d_physics.get_host_ptr()->num_groups() );
        d_db->set("spectral_shape",spectrum);

        // make the uniform source
        auto source = std::make_shared<cuda_mc::Uniform_Source<Geometry> >(
                d_db, d_geometry);
        source->build_source(d_shape);

        // make the solver
        d_fixed_solver = std::make_shared<Fixed_Source_Solver_t>();

        // set it
        d_fixed_solver->set(transporter, source);

        // assign the base solver
        d_solver = d_fixed_solver;
    }
    else
    {
        throw profugus::assertion(
            "Undefined problem type; choose eigenvalue or fixed");
    }

    ENSURE(d_solver);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the problem.
 */
template <class Geometry>
void Manager_Cuda<Geometry>::solve()
{
    if (d_db->template get<bool>("do_transport", true))
    {
        SCOPED_TIMER("Manager_Cuda.solve");

        SCREEN_MSG("Executing solver");

        // run the appropriate solver
#if 0
        if (d_keff_solver)
        {
            CHECK(!d_fixed_solver);
            d_keff_solver->solve();
        }
        else
#endif
        if (d_fixed_solver)
        {
            //CHECK(!d_keff_solver);
            d_fixed_solver->solve();
        }
        else
        {
            throw profugus::assertion("No legitimate solver built");
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do output.
 */
template <class Geometry>
void Manager_Cuda<Geometry>::output()
{
    using std::string;

    SCOPED_TIMER("Manager_Cuda.output");

    SCREEN_MSG("Outputting data");

    // >>> OUTPUT FINAL DATABASE

    // output the final database
    if (d_node == 0)
    {
        std::ostringstream m;
        m << d_problem_name << "_db.xml";
        Teuchos::writeParameterListToXmlFile(*d_db, m.str());
    }

    profugus::global_barrier();

    // >>> OUTPUT SOLUTION (only available if HDF5 is on)
#ifdef USE_HDF5

    // Output filename
    std::ostringstream m;
    m << d_problem_name << "_output.h5";
    string outfile = m.str();

    // scalar output for kcode
#if 0
    if (d_keff_solver)
    {
        // get the kcode tally
        auto keff = d_keff_solver->keff_tally();
        CHECK(keff);

        // make the hdf5 file
        profugus::Serial_HDF5_Writer writer;
        writer.open(outfile);

        // output scalar quantities
        writer.begin_group("keff");

        writer.write(string("mean"), keff->mean());
        writer.write(string("variance"), keff->variance());
        writer.write(string("num_active_cycles"),
                     static_cast<int>(keff->cycle_count()));
        writer.write(string("cycle_estimates"), keff->all_keff());

        writer.end_group();

        writer.close();
    }
#endif

#endif // USE_HDF5
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build geometry
 */
template <class Geometry>
void Manager_Cuda<Geometry>::build_geometry(RCP_ParameterList master_db)
{
    INSIST( master_db->isSublist("MESH"), "Must have MESH sublist." );

    const auto &mesh_db = Teuchos::sublist(master_db, "MESH");

    // Ensure all required parameters are present
    REQUIRE( mesh_db->isType<Array_Dbl>("x_edges") );
    REQUIRE( mesh_db->isType<Array_Dbl>("y_edges") );
    REQUIRE( mesh_db->isType<Array_Dbl>("z_edges") );
    REQUIRE( mesh_db->isType<Array_Int>("matids") );

    const auto &x_edges = mesh_db->get<Array_Dbl>("x_edges");
    const auto &y_edges = mesh_db->get<Array_Dbl>("y_edges");
    const auto &z_edges = mesh_db->get<Array_Dbl>("z_edges");
    const auto &matids  = mesh_db->get<Array_Int>("matids");

    // Build Mesh
    auto geometry_host = std::make_shared<cuda_profugus::Mesh_Geometry>(
        x_edges.toVector(),y_edges.toVector(),z_edges.toVector());

    REQUIRE( matids.size() == ( (x_edges.size()-1) *
                                (y_edges.size()-1) *
                                (z_edges.size()-1) ) );
    geometry_host->set_matids(matids.toVector());

    d_geometry = SDP_Geometry( geometry_host );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build physics
 */
template <class Geometry>
void Manager_Cuda<Geometry>::build_physics(RCP_ParameterList master_db)
{
    typedef profugus::XS_Builder::Matid_Map Matid_Map;
    typedef profugus::XS_Builder::RCP_XS    RCP_XS;

    INSIST(master_db->isSublist("MATERIAL"), "Need MATERIAL database.");

    auto matdb = Teuchos::sublist(master_db,"MATERIAL");

    INSIST(matdb->isType<Array_Str>("mat list"), "Need material list.");
    INSIST(matdb->isType<std::string>("xs library"),
              "Inline cross sections not implemented yet.");

    // get the material list off of the database
    const auto &mat_list = matdb->get<Array_Str>("mat list");

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
    builder.open_and_broadcast(matdb->get<std::string>("xs library"));

    // get the number of groups and moments in the cross section data
    int Ng_data = builder.num_groups();
    int N_data  = builder.pn_order();

    // get the number of groups required
    int g_first = d_db->get("g_first", 0);
    int g_last  = d_db->get("g_last", Ng_data - 1);
    INSIST( (1 + (g_last - g_first) <= Ng_data),
            "Energy group range exceeds number of groups in data." );

    // build the cross sections (always build P0 for Monte Carlo)
    builder.build(matids, 0, g_first, g_last);
    RCP_XS xs = builder.get_xs();
    CHECK(xs->num_mat() == matids.size());
    CHECK(xs->num_groups() == 1 + (g_last - g_first));
    d_db->set("num_groups",xs->num_groups());

    // Cuda physics also needs device pointer
    d_xs_dev = cuda::shared_device_ptr<cuda_profugus::XS_Device>(*xs);

    // make the physics
    auto physics_host = std::make_shared<Physics_t>(d_db, xs, d_xs_dev);

    // set the geometry in the physics
    physics_host->set_geometry(d_geometry);

    // Make shared device pointer
    d_physics = SDP_Physics( physics_host );
}


} // end namespace cuda_mc

#endif // cuda_mc_Manager_Cuda_t_cuh

//---------------------------------------------------------------------------//
//                 end of Manager_Cuda.t.cuh
//---------------------------------------------------------------------------//
