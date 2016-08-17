//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc_driver/Problem_Builder.t.cuh
 * \author Thomas M. Evans
 * \date   Wed Mar 12 22:25:22 2014
 * \brief  Problem_Builder member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_driver_Problem_Builder_t_hh
#define cuda_mc_driver_Problem_Builder_t_hh

#include <utility>
#include <algorithm>
#include <sstream>

#include "comm/global.hh"
#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "utils/String_Functions.hh"
#include "xs/XS_Builder.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "cuda_utils/Hardware.hh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "cuda_mc/Box_Shape.hh"
#include "cuda_mc/Cell_Tally.hh"
#include "cuda_mc/VR_Analog.hh"
#include "cuda_mc/VR_Roulette.hh"
#include "Problem_Builder.hh"
#include "Geometry_Builder.hh"

namespace cuda_mc
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

    // Get number of visible GPUs
    int num_devices;
    auto err = cudaGetDeviceCount(&num_devices);
    CHECK( cudaSuccess == err );

    // Assign MPI tasks to devices by MPI rank
    int device_id = d_db->get("device_id",0);
    cuda_utils::Hardware<cuda_utils::arch::Device>::set_device( device_id );
    CHECK( cuda_utils::Hardware<cuda_utils::arch::Device>::have_acquired() );

    // Set the default CUDA block size.
    if(d_db->isParameter("block_size"))
    {
        cuda_utils::Hardware<cuda_utils::arch::Device>::set_default_block_size(
            d_db->get<int>("block_size") );
    }

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
    INSIST(d_matdb->isParameter("xs library"),
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
    INSIST(1 + (g_last - g_first) <= Ng_data, 
           "Energy group range exceeds number of groups in data" );

    // build the cross sections (always build P0 for Monte Carlo)
    builder.build(matids, 0, g_first, g_last);
    RCP_XS xs = builder.get_xs();
    CHECK(xs->num_mat() == matids.size());
    CHECK(xs->num_groups() == 1 + (g_last - g_first));

    // make the physics
    auto physics = std::make_shared<Physics_t>(*d_db, *xs);

    // set the geometry in the physics
    physics->set_geometry(d_geometry);

    // copy physics to the device.
    d_physics = SDP_Physics( physics );
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
            std::make_shared<cuda_profugus::VR_Roulette<Geom_t> >(*d_db);
    }
    else if (to_lower(var) == "analog")
    {
        d_var_reduction =
            std::make_shared<cuda_profugus::VR_Analog<Geom_t> >();
    }
    else
    {
        INSIST( to_lower(var) == "roulette" ||
                to_lower(var) == "analog",
                "Variance reduction type unknown" );
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
    auto box_shape = std::make_shared<cuda_profugus::Box_Shape>(
        box[0], box[1], box[2], box[3], box[4], box[5]);
    d_shape = SDP_Shape( box_shape );

    // get the source spectrum and add it to the main database
    const auto &shape = source_db.get<OneDArray_dbl>("spectrum");
    CHECK(shape.size() == d_physics.get_host_ptr()->num_groups());

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
    using cuda_profugus::Mesh_Geometry;

    // make the tallier
    d_tallier = std::make_shared<Tallier_t>();

    // add the cell tally.
    if ( d_db->isSublist("cell_tally_db") )
    {
        int num_batch = d_db->get<int>("num_batch",1);
        auto cell_tally = std::make_shared<cuda_profugus::Cell_Tally<Geometry> >(
            d_db, d_geometry, num_batch );
        d_tallier->add_pathlength_tally( cell_tally );
    }

    ENSURE(d_tallier);
}

//---------------------------------------------------------------------------//
/*!
 * \build the SPN problem
 */
template <class Geometry>
void Problem_Builder<Geometry>::build_spn_problem()
{ /* ... */ }

} // end namespace cuda_mc

#endif // cuda_mc_driver_Problem_Builder_t_cuh

//---------------------------------------------------------------------------//
//                 end of Problem_Builder.t.hh
//---------------------------------------------------------------------------//
