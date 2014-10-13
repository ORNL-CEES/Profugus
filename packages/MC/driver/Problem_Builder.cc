//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   driver/Problem_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 12 22:25:22 2014
 * \brief  Problem_Builder member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <utility>
#include <algorithm>
#include <sstream>

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "utils/String_Functions.hh"
#include "xs/XS_Builder.hh"
#include "geometry/Mesh_Geometry.hh"
#include "mc/Box_Shape.hh"
#include "mc/VR_Analog.hh"
#include "mc/VR_Roulette.hh"
#include "mc/Fission_Matrix_Tally.hh"
#include "Problem_Builder.hh"

namespace mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Problem_Builder::Problem_Builder()
    : d_comm(Teuchos::DefaultComm<int>::getComm())
{
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Setup the problem.
 */
void Problem_Builder::setup(const std::string &xml_file)
{
    // make the master parameterlist
    auto master = Teuchos::rcp(new ParameterList(""));

    // read the data on every domain
    Teuchos::updateParametersFromXmlFileAndBroadcast(
        xml_file.c_str(), master.ptr(), *d_comm);

    // validate the parameter list
    INSIST(master->isSublist("CORE"),
            "CORE block not defined in input.");
    INSIST(master->isSublist("ASSEMBLIES"),
            "ASSEMBLIES block not defined in input.");
    INSIST(master->isSublist("PINS"),
            "PINS block not defined in input.");
    INSIST(master->isSublist("MATERIAL"),
            "MATERIAL block not defined in input.");
    INSIST(master->isSublist("PROBLEM"),
            "PROBLEM block not defined in input.");

    // store the individual parameter lists
    d_coredb   = Teuchos::sublist(master, "CORE");
    d_assblydb = Teuchos::sublist(master, "ASSEMBLIES");
    d_pindb    = Teuchos::sublist(master, "PINS");
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
    build_geometry();

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

    // build the tallier
    build_tallies();
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Build the Monte Carlo geometry.
 */
void Problem_Builder::build_geometry()
{
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
    d_geometry = std::make_shared<Geometry_t>(core);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build an axial lattice.
 */
auto Problem_Builder::build_axial_lattice(
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
                                    ids, r, pin_matid, pitch, height));

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

//---------------------------------------------------------------------------//
/*!
 * \brief Build the Monte Carlo MG physics.
 *
 * For now, we build all the cross sections in the problem on every domain.
 * For the mini-app, this is not expected to be an overburdening cost.
 */
void Problem_Builder::build_physics()
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
void Problem_Builder::build_var_reduction()
{
    using profugus::to_lower;

    REQUIRE(!d_db.is_null());

    // the default is to do roulette
    const auto &var = to_lower(
        d_db->get<std::string>("variance reduction", std::string("roulette")));

    // build the appropriate variance reduction
    if (to_lower(var) == "roulette")
    {
        d_var_reduction = std::make_shared<profugus::VR_Roulette>(d_db);
    }
    else if (to_lower(var) == "analog")
    {
        d_var_reduction = std::make_shared<profugus::VR_Analog>();
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
void Problem_Builder::build_source(const ParameterList &source_db)
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
void Problem_Builder::build_tallies()
{
    using profugus::Mesh_Geometry;
    using profugus::Fission_Matrix_Tally;

    // make the tallier
    d_tallier = std::make_shared<Tallier_t>();
    d_tallier->set(d_geometry, d_physics);

    // check for fission matrix tallies
    if (d_db->isSublist("fission_matrix_db"))
    {
        // get the database
        const ParameterList &fdb = d_db->sublist("fission_matrix_db");

        // validate that there is a fission matrix mesh defined
        VALIDATE(fdb.isParameter("x_bounds"), "Failed to define x-boundaries "
                 << "for fission matrix mesh");
        VALIDATE(fdb.isParameter("y_bounds"), "Failed to define y-boundaries "
                 << "for fission matrix mesh");
        VALIDATE(fdb.isParameter("z_bounds"), "Failed to define z-boundaries "
                 << "for fission matrix mesh");

        // get the fission matrix mesh boundaries
        auto xb = fdb.get<OneDArray_dbl>("x_bounds").toVector();
        auto yb = fdb.get<OneDArray_dbl>("y_bounds").toVector();
        auto zb = fdb.get<OneDArray_dbl>("z_bounds").toVector();
        Fission_Matrix_Tally::SP_Mesh_Geometry geo(
            std::make_shared<Mesh_Geometry>(xb, yb, zb));

        // build the fission matrix
        auto fm_tally(std::make_shared<Fission_Matrix_Tally>(
                          d_db, d_physics, geo));
        CHECK(fm_tally);

        // add this to the tallier
        d_tallier->add_pathlength_tally(fm_tally);
        d_tallier->add_source_tally(fm_tally);
    }

    ENSURE(d_tallier);
}

} // end namespace mc

//---------------------------------------------------------------------------//
//                 end of Problem_Builder.cc
//---------------------------------------------------------------------------//
