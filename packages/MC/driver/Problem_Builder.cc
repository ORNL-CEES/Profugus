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

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "xs/XS_Builder.hh"
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
    Insist (master->isSublist("CORE"),
            "CORE block not defined in input.");
    Insist (master->isSublist("ASSEMBLIES"),
            "ASSEMBLIES block not defined in input.");
    Insist (master->isSublist("PINS"),
            "PINS block not defined in input.");
    Insist (master->isSublist("MATERIAL"),
            "MATERIAL block not defined in input.");
    Insist (master->isSublist("PROBLEM"),
            "PROBLEM block not defined in input.");

    // store the individual parameter lists
    d_coredb   = Teuchos::sublist(master, "CORE");
    d_assblydb = Teuchos::sublist(master, "ASSEMBLIES");
    d_pindb    = Teuchos::sublist(master, "PINS");
    d_matdb    = Teuchos::sublist(master, "MATERIAL");
    d_db       = Teuchos::sublist(master, "PROBLEM");

    Check (!d_db.is_null());

    // build the geometry
    build_geometry();

    // build build physics
    build_physics();

    // build the external source (there won't be one for k-eigenvalue
    // problems)
    if (master->isSublist("SOURCE"))
    {
        build_source(master->sublist("SOURCE"));
    }
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Build the Monte Carlo geometry.
 */
void Problem_Builder::build_geometry()
{
    Require (d_coredb->isParameter("axial list"));
    Require (d_coredb->isParameter("axial height"));
    Require (d_assblydb->isParameter("assembly list"));

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
    Check (axial_height.size() == num_axial);

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
        Check (d_coredb->isParameter(axial_list[k]));

        // get the core map at this axial level
        const auto &core_map = d_coredb->get<TwoDArray_int>(axial_list[k]);
        Check (core_map.getNumCols() == num_x);
        Check (core_map.getNumRows() == num_y);

        // run through the core-map on this level and build all of the
        // assemblies
        for (int j = 0; j < num_y; ++j)
        {
            for (int i = 0; i < num_x; ++i)
            {
                // get the id of this assembly (in range [0, N))
                aid = core_map(j, i);
                Check (aid < assbly_list.size());
                Check (d_assblydb->isParameter(assbly_list[aid]));

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

                Check (assemblies[k].count(aid));
            }
        }
        Check (assemblies[k].size() == a2rtk[k].size());

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
            Check (a.second);
            Check (a2rtk[k][a.first] + axial_offset >= 0 &&
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
                Check (a2rtk[k][map(j, i)] + axial_offset >= 0 &&
                       a2rtk[k][map(j, i)] + axial_offset < num_assemblies);
                core->id(i, j, k) = a2rtk[k][map(j, i)] + axial_offset;
            }
        }

        // update the axial offset
        axial_offset += assemblies[k].size();
    }
    Check (axial_offset == num_assemblies);

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
    Require (d_pindb->isParameter("pin list"));

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
            Check (d_pindb->isSublist(pin_list[pid]));

            // build the pin if we haven't already
            if (!pins.count(pid))
            {
                // get the pin-sublist
                const auto &pindb = d_pindb->sublist(pin_list[pid]);
                Check (pindb.isParameter("pitch"));
                Check (pindb.isParameter("matid"));

                // get the pin pitch and overall pin matid
                double pitch  = pindb.get<double>("pitch");
                int pin_matid = pindb.get<int>("matid");

                // make empty inner cylinders
                Pin_Cell_t::Vec_Dbl r;
                Pin_Cell_t::Vec_Int ids;

                // see if the pin has any internal cylinders
                if (pindb.isParameter("radii"))
                {
                    Check (pindb.isParameter("radial matids"));
                    const auto &radii  = pindb.get<OneDArray_dbl>("radii");
                    const auto &matids = pindb.get<OneDArray_int>(
                        "radial matids");

                    r.insert(r.end(), radii.begin(), radii.end());
                    ids.insert(ids.end(), matids.begin(), matids.end());

                    Check (r.size() == radii.size());
                    Check (ids.size() == matids.size());
                }

                // build the pin
                SP_Pin_Cell pin(std::make_shared<Pin_Cell_t>(
                                    ids, r, pin_matid, pitch, 1.0));

                // add it
                pins.emplace(pid, pin);

                // make the input-to-rtk pin id
                inp2rtk_id.emplace(pid, pins.size() - 1);
            }
            Check (pins.count(pid));
        }
    }
    Check (inp2rtk_id.size() == pins.size());

    // make the lattice geometry
    SP_Lattice lattice(std::make_shared<Lattice_t>(
                           num_x, num_y, 1, pins.size()));

    // assign the pins to the lattice
    for (const auto &p : pins)
    {
        Check (p.second);
        Check (inp2rtk_id[p.first] >= 0 && inp2rtk_id[p.first] < pins.size());
        lattice->assign_object(p.second, inp2rtk_id[p.first]);
    }

    // assign the RTK pin ids to the lattice
    for (int j = 0; j < num_y; ++j)
    {
        for (int i = 0; i < num_x; ++i)
        {
            Check (inp2rtk_id[map(j, i)] >= 0 &&
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

    Require (d_matdb->isParameter("mat list"));
    Validate (d_matdb->isParameter("xs library"),
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
    Check (matids.size() == mat_list.size());

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
    Validate (1 + (g_last - g_first) <= Ng_data, "Energy group range exceeds "
              << "number of groups in data, 1 + g_last - g_first = "
              << 1 + (g_last - g_first) << " > " << Ng_data);

    // build the cross sections (always build P0 for Monte Carlo)
    builder.build(matids, 0, g_first, g_last);
    RCP_XS xs = builder.get_xs();
    Check (xs->num_mat() == matids.size());
    Check (xs->num_groups() == 1 + (g_last - g_first));

    // make the physics
    d_physics = std::make_shared<Physics_t>(d_db, xs);

    // set the geometry in the physics
    d_physics->set_geometry(d_geometry);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the external source.
 *
 *\param source_db
 */
void Problem_Builder::build_source(const ParameterList &source_db)
{
}

} // end namespace mc

//---------------------------------------------------------------------------//
//                 end of Problem_Builder.cc
//---------------------------------------------------------------------------//
