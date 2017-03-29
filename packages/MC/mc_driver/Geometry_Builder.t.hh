//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_driver/Geometry_Builder.t.hh
 * \author Steven Hamilton
 * \date   Wed Nov 25 12:58:58 2015
 * \brief  Geometry_Builder template member definitions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_driver_Geometry_Builder_t_hh
#define mc_driver_Geometry_Builder_t_hh

#include "Geometry_Builder.hh"

namespace mc
{

auto Geometry_Builder<profugus::Core>::build(
    RCP_ParameterList master) -> SP_Geometry
{
    // validate the parameter list
    INSIST(master->isSublist("PROBLEM"),
            "PROBLEM block not defined in input.");
    INSIST(master->isSublist("CORE"),
            "CORE block not defined in input.");
    INSIST(master->isSublist("ASSEMBLIES"),
            "ASSEMBLIES block not defined in input.");
    INSIST(master->isSublist("PINS"),
            "PINS block not defined in input.");

    // store the individual parameter lists
    d_db       = Teuchos::sublist(master, "PROBLEM");
    d_coredb   = Teuchos::sublist(master, "CORE");
    d_assblydb = Teuchos::sublist(master, "ASSEMBLIES");
    d_pindb    = Teuchos::sublist(master, "PINS");

    REQUIRE(d_coredb->isParameter("axial list"));
    REQUIRE(d_coredb->isParameter("axial height"));
    REQUIRE(d_assblydb->isParameter("assembly list"));

    // get the axial list and heights
    const auto &axial_list   = d_coredb->get<OneDArray_str>("axial list");
    const auto &axial_height = d_coredb->get<OneDArray_dbl>("axial height");

    // build the core (all axial core maps have the same radial dimensions, so
    // we can just use the first core map here) -> calculate the maximum
    // number of objects we may have (most of these will probably not be
    // needed)
    const auto &base_map = d_coredb->get<TwoDArray_int>(axial_list[0]);

    // get the core dimensions (radially in assemblies, axially in levels);
    // remember the twoD arrays are entered [j][i] (i moves fastest in
    // COLUMN-MAJOR---FORTRAN---style, so it goes in the column index)
    int num_x     = base_map.getNumCols();
    int num_y     = base_map.getNumRows();
    int num_axial = axial_list.size();
    CHECK(axial_height.size() == num_axial);

    // get the assembly list
    const auto &assbly_list = d_assblydb->get<OneDArray_str>("assembly list");

    // make a hash of the assemblies that need to be built; we build unique
    // assemblies on each *axial* level
    std::vector<Lattice_Hash> assemblies(num_axial);
    int                       aid = 0;

    // map of assembly ids to rtk-ids on each axial level
    std::vector<std::unordered_map<int, int>> a2rtk(num_axial);

    // number of unique assemblies in core
    int num_assemblies = 0;

    // iterate through the core and build each assembly by axial level
    for (int k = 0; k < num_axial; ++k)
    {
        CHECK(d_coredb->isParameter(axial_list[k]));

        // get the core map at this axial level
        const auto &core_map = d_coredb->get<TwoDArray_int>(axial_list[k]);
        CHECK(core_map.getNumCols() == num_x);
        CHECK(core_map.getNumRows() == num_y);

        // run through the core-map on this level and build all of the
        // assemblies
        for (int j = 0; j < num_y; ++j)
        {
            for (int i = 0; i < num_x; ++i)
            {
                // get the id of this assembly (in range [0, N))
                aid = core_map(j, i);
                CHECK(aid < assbly_list.size());
                CHECK(d_assblydb->isParameter(assbly_list[aid]));

                // get the assembly map for this lattice
                const auto &assbly_map = d_assblydb->get<TwoDArray_int>(
                    assbly_list[aid]);

                // build the lattice if we have not added it already
                if (!assemblies[k].count(aid))
                {
                    SP_Lattice lattice = build_axial_lattice(
                        assbly_map, axial_height[k]);

                    // add it to the hash table
                    assemblies[k].emplace(aid, lattice);

                    // map the assembly id to rtk id
                    a2rtk[k].emplace(aid, assemblies[k].size() - 1);
                }

                CHECK(assemblies[k].count(aid));
            }
        }
        CHECK(assemblies[k].size() == a2rtk[k].size());

        // add up the number of unique assemblies (on each axial level)
        num_assemblies += assemblies[k].size();
    }

    // make the core
    SP_Core core(std::make_shared<Core_t>(
                     num_x, num_y, num_axial, num_assemblies));

    // now add the assemblies to the core
    int axial_offset = 0;
    for (int k = 0; k < num_axial; ++k)
    {
        // assign the assemblies
        for (const auto &a : assemblies[k])
        {
            CHECK(a.second);
            CHECK(a2rtk[k][a.first] + axial_offset >= 0 &&
                   a2rtk[k][a.first] + axial_offset < num_assemblies);
            core->assign_object(a.second, a2rtk[k][a.first] + axial_offset);
        }

        // get the core map at this axial level
        const auto &map = d_coredb->get<TwoDArray_int>(axial_list[k]);

        // assign the RTK pin ids to the core
        for (int j = 0; j < num_y; ++j)
        {
            for (int i = 0; i < num_x; ++i)
            {
                CHECK(a2rtk[k][map(j, i)] + axial_offset >= 0 &&
                       a2rtk[k][map(j, i)] + axial_offset < num_assemblies);
                core->id(i, j, k) = a2rtk[k][map(j, i)] + axial_offset;
            }
        }

        // update the axial offset
        axial_offset += assemblies[k].size();
    }
    CHECK(axial_offset == num_assemblies);

    // set the boundary conditions
    def::Vec_Int boundary(6, 0);
    if (d_db->get<std::string>("boundary") == "reflect")
    {
        CHECK(d_db->isSublist("boundary_db"));
        CHECK(d_db->sublist("boundary_db").isParameter("reflect"));
        const auto &bnd_array =
            d_db->sublist("boundary_db").get<OneDArray_int>("reflect");
        std::copy(bnd_array.begin(), bnd_array.end(), boundary.begin());
    }
    core->set_reflecting(boundary);

    // complete the core
    core->complete(0.0, 0.0, 0.0);

    // make the geometry
    return std::make_shared<Geom_t>(core);
}

auto Geometry_Builder<profugus::Core>::build_axial_lattice(
        const TwoDArray_int &map,
        double               height) -> SP_Lattice
{
    REQUIRE(d_pindb->isParameter("pin list"));

    // get the pin list from the pindb
    const auto &pin_list = d_pindb->get<OneDArray_str>("pin list");

    // make a hash of unique pins to build in this assembly
    Pin_Hash pins;
    int      pid = 0;

    // number of pins (remember the indices are flipped because the data is
    // loaded COLUMN-MAJOR)
    int num_x = map.getNumCols();
    int num_y = map.getNumRows();

    // map of input-pin-id to local (RTK) pin-id; the RTK pin-id must be in
    // the range [0,num_pins) whereas the input-pin-id can be any index
    std::unordered_map<int, int> inp2rtk_id;

    // build the pins
    for (int j = 0; j < num_y; ++j)
    {
        for (int i = 0; i < num_x; ++i)
        {
            // get the pinid
            pid = map(j, i);
            CHECK(d_pindb->isSublist(pin_list[pid]));

            // build the pin if we haven't already
            if (!pins.count(pid))
            {
                // get the pin-sublist
                const auto &pindb = d_pindb->sublist(pin_list[pid]);
                CHECK(pindb.isParameter("pitch"));
                CHECK(pindb.isParameter("matid"));

                // get the pin pitch and overall pin matid
                double pitch  = pindb.get<double>("pitch");
                int pin_matid = pindb.get<int>("matid");
                int num_segments = 1;
                if (pindb.isType<int>("num segments"))
                    num_segments = pindb.get<int>("num segments");
                CHECK(num_segments == 1 || num_segments == 4);

                // make empty inner cylinders
                Pin_Cell_t::Vec_Dbl r;
                Pin_Cell_t::Vec_Int ids;

                // see if the pin has any internal cylinders
                if (pindb.isParameter("radii"))
                {
                    CHECK(pindb.isParameter("radial matids"));
                    const auto &radii  = pindb.get<OneDArray_dbl>("radii");
                    const auto &matids = pindb.get<OneDArray_int>(
                        "radial matids");

                    r.insert(r.end(), radii.begin(), radii.end());
                    ids.insert(ids.end(), matids.begin(), matids.end());

                    CHECK(r.size() == radii.size());
                    CHECK(ids.size() == matids.size());
                }

                // build the pin
                SP_Pin_Cell pin(std::make_shared<Pin_Cell_t>(
                    ids, r, pin_matid, pitch, height, num_segments));

                // add it
                pins.emplace(pid, pin);

                // make the input-to-rtk pin id
                inp2rtk_id.emplace(pid, pins.size() - 1);
            }
            CHECK(pins.count(pid));
        }
    }
    CHECK(inp2rtk_id.size() == pins.size());

    // make the lattice geometry
    SP_Lattice lattice(std::make_shared<Lattice_t>(
                           num_x, num_y, 1, pins.size()));

    // assign the pins to the lattice
    for (const auto &p : pins)
    {
        CHECK(p.second);
        CHECK(inp2rtk_id[p.first] >= 0 && inp2rtk_id[p.first] < pins.size());
        lattice->assign_object(p.second, inp2rtk_id[p.first]);
    }

    // assign the RTK pin ids to the lattice
    for (int j = 0; j < num_y; ++j)
    {
        for (int i = 0; i < num_x; ++i)
        {
            CHECK(inp2rtk_id[map(j, i)] >= 0 &&
                   inp2rtk_id[map(j, i)] < pins.size());
            lattice->id(i, j, 0) = inp2rtk_id[map(j, i)];
        }
    }

    lattice->complete(0.0, 0.0, 0.0);
    return lattice;
}

auto Geometry_Builder<profugus::Mesh_Geometry>::build(
    RCP_ParameterList master) -> SP_Geometry
{
    auto problem_db = Teuchos::sublist(master, "PROBLEM");
    auto mesh_db    = Teuchos::sublist(master, "MESH");

    // Ensure all required parameters are present
    REQUIRE( mesh_db->isParameter("x_edges") );
    REQUIRE( mesh_db->isParameter("y_edges") );
    REQUIRE( mesh_db->isParameter("z_edges") );
    REQUIRE( mesh_db->isParameter("matids") );

    auto x_edges = mesh_db->get<OneDArray_dbl>("x_edges");
    auto y_edges = mesh_db->get<OneDArray_dbl>("y_edges");
    auto z_edges = mesh_db->get<OneDArray_dbl>("z_edges");
    auto matids  = mesh_db->get<OneDArray_int>("matids");

    // Build Mesh
    auto geom = std::make_shared<profugus::Mesh_Geometry>(
        x_edges.toVector(),y_edges.toVector(),z_edges.toVector());

    // Convert matids to SP<Vec_Int> and pass to geometry
    auto sp_matids = std::make_shared<def::Vec_Int>(
        matids.begin(),matids.end());

    REQUIRE( sp_matids->size() == ( (x_edges.size()-1) *
                                    (y_edges.size()-1) *
                                    (z_edges.size()-1) ) );
    geom->set_matids(sp_matids);

    // set the boundary conditions
    def::Vec_Int boundary(6, 0);
    if (problem_db->get<std::string>("boundary") == "reflect")
    {
        CHECK(problem_db->isSublist("boundary_db"));
        CHECK(problem_db->sublist("boundary_db").isParameter("reflect"));
        const auto &bnd_array =
            problem_db->sublist("boundary_db").get<OneDArray_int>("reflect");
        std::copy(bnd_array.begin(), bnd_array.end(), boundary.begin());
    }
    geom->set_reflecting(boundary);

    return geom;
}

} // end namespace mc

#endif // mc_driver_Geometry_Builder_t_hh

//---------------------------------------------------------------------------//
//                 end of Geometry_Builder.t.hh
//---------------------------------------------------------------------------//
