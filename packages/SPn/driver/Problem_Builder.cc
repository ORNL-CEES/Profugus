//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   driver/Problem_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 12 22:25:22 2014
 * \brief  Problem_Builder member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>
#include <numeric>
#include <utility>

#include "Teuchos_XMLParameterListHelpers.hpp"
#ifdef COMM_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#else
#include "Teuchos_DefaultSerialComm.hpp"
#endif

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "xs/XS_Builder.hh"
#include "mesh/Partitioner.hh"
#include "Problem_Builder.hh"
#include "comm/global.hh"

namespace spn
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Problem_Builder::Problem_Builder()
{
#ifdef COMM_MPI
    d_comm = Teuchos::rcp( new Teuchos::MpiComm<int>(profugus::communicator) );
#else
    d_comm = Teuchos::rcp( new Teuchos::SerialComm<int>() );
#endif
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
    INSIST(master->isSublist("MATERIAL"),
            "MATERIAL block not defined in input.");
    INSIST(master->isSublist("PROBLEM"),
            "PROBLEM block not defined in input.");

    // store the individual parameter lists
    d_coredb   = Teuchos::sublist(master, "CORE");
    d_assblydb = Teuchos::sublist(master, "ASSEMBLIES");
    d_matdb    = Teuchos::sublist(master, "MATERIAL");
    d_db       = Teuchos::sublist(master, "PROBLEM");

    // kill any remaining sources
    d_source = RCP_Source();
    CHECK(d_source.is_null());

    CHECK(!d_db.is_null());

    // build mesh
    build_mesh();

    // build the material ids on the mesh
    build_matids();

    // build material database
    build_matdb();

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
 * \brief Build/partition the mesh.
 */
void Problem_Builder::build_mesh()
{
    using def::I; using def::J; using def::K;

    REQUIRE(d_coredb->isParameter("axial list"));
    REQUIRE(d_coredb->isParameter("axial height"));
    REQUIRE(d_assblydb->isParameter("assembly list"));
    REQUIRE(d_assblydb->isParameter("pin pitch"));
    REQUIRE(d_db->isParameter("radial mesh"));
    REQUIRE(d_db->isParameter("axial mesh"));
    REQUIRE(d_db->isParameter("symmetry"));

    // get the axial core map and heights
    const auto &axial_list   = d_coredb->get<OneDArray_str>("axial list");
    const auto &axial_height = d_coredb->get<OneDArray_dbl>("axial height");
    CHECK(!axial_list.empty());
    CHECK(axial_list.size() == axial_height.size());

    // build the mesh dimensions (all axial core maps have the same radial
    // dimensions, so we can just use the first core map here)
    const auto &core_map = d_coredb->get<TwoDArray_int>(axial_list[0]);

    // get the core dimensions (radially in assemblies, axially in levels);
    // remember the twoD arrays are entered [j][i] (i moves fastest in
    // COLUMN-MAJOR---FORTRAN---style, so it goes in the column index)
    int num_axial_levels = axial_list.size();
    d_Na[I]              = core_map.getNumCols();
    d_Na[J]              = core_map.getNumRows();

    // all assemblies have the same radial dimensions, so use the first one to
    // get the core dimensions
    const auto &assbly_list = d_assblydb->get<OneDArray_str>("assembly list");
    const auto &assbly_map  = d_assblydb->get<TwoDArray_int>(assbly_list[0]);

    // get pins (same caveats on ordering as for the core map)
    d_Np[I] = assbly_map.getNumCols();
    d_Np[J] = assbly_map.getNumRows();

    // get the pin pitch
    double pitch = d_assblydb->get<double>("pin pitch");
    CHECK(pitch > 0.0);

    // get the core dimensions
    double dx = d_Na[I] * d_Np[I] * pitch;
    double dy = d_Na[J] * d_Np[J] * pitch;
    double dz = std::accumulate(axial_height.begin(), axial_height.end(), 0.0);

    // get the mesh dimensions
    int radial_mesh        = d_db->get<int>("radial mesh");
    const auto &axial_mesh = d_db->get<OneDArray_int>("axial mesh");
    CHECK(axial_mesh.size() == axial_height.size());

    // set the mesh radial mesh dimensions
    int ncx = radial_mesh * d_Np[I] * d_Na[I];
    int ncy = radial_mesh * d_Np[J] * d_Na[J];
    d_db->set("num_cells_i", ncx);
    d_db->set("delta_x", pitch);
    d_db->set("num_cells_j", ncy);
    d_db->set("delta_y", pitch);

    // set the mesh axial dimensions
    int ncz = std::accumulate(axial_mesh.begin(), axial_mesh.end(), 0);
    OneDArray_dbl z_edges(ncz + 1, 0.0);

    // iterate through axial levels and build the mesh
    int ctr      = 0;
    double delta = 0.0;
    for (int k = 0; k < num_axial_levels; ++k)
    {
        // get the width of cells in this level
        delta = axial_height[k] / axial_mesh[k];

        // iterate through axial cells on this level
        for (int n = 0; n < axial_mesh[k]; ++n)
        {
            CHECK(ctr + 1 < z_edges.size());
            z_edges[ctr+1] = z_edges[ctr] + delta;
            ++ctr;
        }
    }
    CHECK(ctr == ncz);
    CHECK(profugus::soft_equiv(dz, z_edges.back(), 1.0e-12));

    // set the axial mesh
    d_db->set("z_edges", z_edges);

    // partition the mesh
    profugus::Partitioner p(d_db);
    p.build();

    // assign mesh objects
    d_mesh    = p.get_mesh();
    d_indexer = p.get_indexer();
    d_gdata   = p.get_global_data();

    ENSURE(!d_mesh.is_null());
    ENSURE(!d_indexer.is_null());
    ENSURE(!d_gdata.is_null());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the material ids.
 */
void Problem_Builder::build_matids()
{
    using def::I; using def::J; using def::K;

    // size the matids array
    d_matids.resize(d_mesh->num_cells());

    // k-mesh levels for each axial level
    int k_begin = 0, k_end = 0;

    // global radial core map
    TwoDArray_int axial_matids(d_gdata->num_cells(J), d_gdata->num_cells(I), 0);

    // get number of cells per axial level
    const auto &axial_mesh = d_db->get<OneDArray_int>("axial mesh");

    // process the axial levels one at a time
    for (int level = 0; level < axial_mesh.size(); ++level)
    {
        // calculate the global matids in this level
        calc_axial_matids(level, axial_matids);

        // determine the begin/end of the axial mesh for this axial level
        k_begin = k_end;
        k_end   = k_begin + axial_mesh[level];
        CHECK(k_end - k_begin == axial_mesh[level]);

        // loop over local cells
        for (int k = k_begin; k < k_end; ++k)
        {
            for (int j = 0; j < d_mesh->num_cells_dim(J); ++j)
            {
                for (int i = 0; i < d_mesh->num_cells_dim(I); ++i)
                {
                    // get the global IJ indices
                    auto global = d_indexer->convert_to_global(i, j);
                    CHECK(global[I] < axial_matids.getNumCols());
                    CHECK(global[J] < axial_matids.getNumRows());

                    // assign the local matid
                    d_matids[d_indexer->l2l(i, j, k)] =
                        axial_matids(global[J], global[I]);
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate material ids in an axial level.
 */
void Problem_Builder::calc_axial_matids(int            level,
                                        TwoDArray_int &matids)
{
    using def::I; using def::J;

    // (x, y) global offsets into the mesh by assembly and pin
    int aoff_x = 0, poff_x = 0;
    int aoff_y = 0, poff_y = 0;

    // get the list of core maps and assembly types
    const auto &axial_list  = d_coredb->get<OneDArray_str>("axial list");
    const auto &assbly_list = d_assblydb->get<OneDArray_str>("assembly list");

    // get the core-map for this axial level
    const auto &core_map = d_coredb->get<TwoDArray_int>(axial_list[level]);
    CHECK(core_map.getNumCols() == d_Na[I]);
    CHECK(core_map.getNumRows() == d_Na[J]);

    // mesh cells per pin
    int radial_mesh = d_db->get<int>("radial mesh");
    CHECK(matids.getNumCols() == d_Na[I] * d_Np[I] * radial_mesh);
    CHECK(matids.getNumRows() == d_Na[J] * d_Np[J] * radial_mesh);

    // loop over all assemblies, get the pin-maps, and assign the material ids
    // to the matids array (remember, all "core arrays" are ordered
    // COLUMN-MAJOR, which means matids[j, i])

    // material id of a pin
    int matid = 0;

    // loop over assemblies in J
    for (int aj = 0; aj < d_Na[J]; ++aj)
    {
        // set the y-offset for this assembly
        aoff_y = (radial_mesh * d_Np[J]) * aj;

        // loop over assemblies in I
        for (int ai = 0; ai < d_Na[I]; ++ai)
        {
            CHECK(core_map(aj, ai) < assbly_list.size());
            CHECK(d_assblydb->isParameter(assbly_list[core_map(aj, ai)]));

            // get the pin-map for this assembly
            const auto &assbly_map = d_assblydb->get<TwoDArray_int>(
                assbly_list[core_map(aj, ai)]);
            CHECK(assbly_map.getNumCols() == d_Np[I]);
            CHECK(assbly_map.getNumRows() == d_Np[J]);

            // set the x-offset for this assembly
            aoff_x = (radial_mesh * d_Np[I]) * ai;

            // loop over pins in J
            for (int pj = 0; pj < d_Np[J]; ++pj)
            {
                // set the y-offset for this pin
                poff_y = aoff_y + radial_mesh * pj;

                for (int pi = 0; pi < d_Np[I]; ++pi)
                {
                    // set the x-offset for this pin
                    poff_x = aoff_x + radial_mesh * pi;

                    // get the material id for this pin
                    matid = assbly_map(pj, pi);
                    CHECK(matid <
                           d_matdb->get<OneDArray_str>("mat list").size());

                    // loop over the mesh cells in this pin
                    for (int j = 0; j < radial_mesh; ++j)
                    {
                        for (int i = 0; i < radial_mesh; ++i)
                        {
                            CHECK(i + poff_x < matids.getNumCols());
                            CHECK(j + poff_y < matids.getNumRows());
                            matids(j + poff_y, i + poff_x) = matid;
                        }
                    }
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the material data.
 *
 * For now, we build all the cross sections in the problem on every domain.
 * For the mini-app, this is not expected to be an overburdening cost.
 */
void Problem_Builder::build_matdb()
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

    // determine the moment order of the problem
    int pn_order = d_db->get("Pn_order", N_data);
    VALIDATE(pn_order <= N_data, "Requested Pn scattering order of "
              << pn_order << " is greater than available data Pn order of "
              << N_data);

    // get the number of groups required
    int g_first = d_db->get("g_first", 0);
    int g_last  = d_db->get("g_last", Ng_data - 1);
    VALIDATE(1 + (g_last - g_first) <= Ng_data, "Energy group range exceeds "
              << "number of groups in data, 1 + g_last - g_first = "
              << 1 + (g_last - g_first) << " > " << Ng_data);

    // build the cross sections
    builder.build(matids, pn_order, g_first, g_last);
    RCP_XS xs = builder.get_xs();
    CHECK(xs->num_mat() == matids.size());
    CHECK(xs->num_groups() == 1 + (g_last - g_first));

    // build the material database
    d_mat = Teuchos::rcp(new profugus::Mat_DB);

    // set the cross sections
    d_mat->set(xs, d_mesh->num_cells());

    // set the matids in the material database
    d_mat->assign(d_matids);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the external source.
 */
void Problem_Builder::build_source(const ParameterList &source_db)
{
    using def::I; using def::J;

    REQUIRE(source_db.isParameter("source list"));
    REQUIRE(source_db.isParameter("source map"));
    REQUIRE(source_db.isParameter("axial source"));
    REQUIRE(!d_mesh.is_null());
    REQUIRE(!d_mat.is_null());

    // get the list of sources
    const auto &source_list  = source_db.get<OneDArray_str>("source list");

    // global radial source map
    Teuchos::TwoDArray<std::pair<int, double>> gsrc(
        d_gdata->num_cells(J), d_gdata->num_cells(I),
        std::pair<int, double>(0, 0.0));

    // get the source-map
    const auto &src_map = source_db.get<TwoDArray_int>("source map");
    CHECK(src_map.getNumCols() == d_Na[I]);
    CHECK(src_map.getNumRows() == d_Na[J]);

    // (x, y) global offsets into the mesh by assembly and pin
    int aoff_x = 0, poff_x = 0;
    int aoff_y = 0, poff_y = 0;

    // mesh cells per pin
    int radial_mesh = d_db->get<int>("radial mesh");

    // source strength and id
    int    src_id       = 0;
    double src_strength = 0.0;

    // loop over assemblies in J
    for (int aj = 0; aj < d_Na[J]; ++aj)
    {
        // set the y-offset for this assembly
        aoff_y = (radial_mesh * d_Np[J]) * aj;

        // loop over assemblies in I
        for (int ai = 0; ai < d_Na[I]; ++ai)
        {
            // store the id of this source
            src_id = src_map(aj, ai);

            // if the src_id >= 0 then there is a source defined here
            if (src_id >= 0)
            {
                CHECK(src_map(aj, ai) < source_list.size());
                CHECK(source_db.isSublist(source_list[src_map(aj, ai)]));

                // get the source map for this assembly
                const auto &a_src_map = source_db.sublist(
                    source_list[src_id]).get<TwoDArray_dbl>("strength");
                CHECK(a_src_map.getNumCols() == d_Np[I]);
                CHECK(a_src_map.getNumRows() == d_Np[J]);

                // set the x-offset for this assembly
                aoff_x = (radial_mesh * d_Np[I]) * ai;

                // loop over pins in J
                for (int pj = 0; pj < d_Np[J]; ++pj)
                {
                    // set the y-offset for this pin
                    poff_y = aoff_y + radial_mesh * pj;

                    for (int pi = 0; pi < d_Np[I]; ++pi)
                    {
                        // set the x-offset for this pin
                        poff_x = aoff_x + radial_mesh * pi;

                        // get the source strength for this pin
                        src_strength = a_src_map(pj, pi);
                        CHECK(src_strength >= 0.0);

                        // loop over the mesh cells in this pin
                        for (int j = 0; j < radial_mesh; ++j)
                        {
                            for (int i = 0; i < radial_mesh; ++i)
                            {
                                CHECK(i + poff_x < gsrc.getNumCols());
                                CHECK(j + poff_y < gsrc.getNumRows());

                                // store the source id (shape id) for this cell
                                gsrc(j + poff_y, i + poff_x).first = src_id;

                                // store the source strength for this cell
                                gsrc(j + poff_y, i + poff_x).second =
                                    src_strength;
                            }
                        }
                    }
                }
            } // process the source
        } // assembly-I
    } // assembly-J

    // build the source shapes
    profugus::Isotropic_Source::Source_Shapes shapes(
        source_list.size(), profugus::Isotropic_Source::Shape(
            d_mat->xs().num_groups(), 0.0));

    // get the first and last groups for the run
    int g_first = d_db->get<int>("g_first");
    int g_last  = d_db->get<int>("g_last");
    CHECK(1 + g_last - g_first == d_mat->xs().num_groups());

    // loop over sources
    int ctr = 0;
    for (auto itr = source_list.begin(); itr != source_list.end(); ++itr, ++ctr)
    {
        CHECK(source_db.sublist(*itr).isParameter("shape"));

        // get the source shapes
        const auto &shape = source_db.sublist(*itr).get<OneDArray_dbl>(
            "shape");
        CHECK(shape.size() <= d_mat->xs().num_groups());

        // loop from first-to last groups and add the shapes
        for (int g = g_first; g <= g_last; ++g)
        {
            CHECK(g < shape.size());
            CHECK(g-g_first < shapes[ctr].size());
            shapes[ctr][g-g_first] = shape[g];
        }
    }
    CHECK(ctr == source_list.size());

    // build the source
    d_source = Teuchos::rcp(new profugus::Isotropic_Source(
                                d_mesh->num_cells()));

    // make the field of source strengths
    profugus::Isotropic_Source::Source_Field strengths(
        d_mesh->num_cells(), 0.0);

    // make the source ids field
    profugus::Isotropic_Source::ID_Field srcids(d_mesh->num_cells(), 0);

    // get the axial source strength and axial mesh
    const auto &axial_mesh = d_db->get<OneDArray_int>("axial mesh");
    const auto &axial_src  = source_db.get<OneDArray_dbl>("axial source");
    CHECK(axial_src.size() == axial_mesh.size());

    // k-mesh levels for each axial level
    int k_begin = 0, k_end = 0;

    // loop over axial levels
    for (int level = 0; level < axial_mesh.size(); ++level)
    {
        // determine the begin/end of the axial mesh for this axial level
        k_begin = k_end;
        k_end   = k_begin + axial_mesh[level];
        CHECK(k_end - k_begin == axial_mesh[level]);

        // loop over local cells
        for (int k = k_begin; k < k_end; ++k)
        {
            for (int j = 0; j < d_mesh->num_cells_dim(J); ++j)
            {
                for (int i = 0; i < d_mesh->num_cells_dim(I); ++i)
                {
                    // get the global IJ indices
                    auto global = d_indexer->convert_to_global(i, j);
                    CHECK(global[I] < gsrc.getNumCols());
                    CHECK(global[J] < gsrc.getNumRows());

                    // assign the source id
                    srcids[d_indexer->l2l(i, j, k)] =
                        gsrc(global[J], global[I]).first;

                    // assign the source strength
                    strengths[d_indexer->l2l(i, j, k)] =
                        gsrc(global[J], global[I]).second * axial_src[level];
                }
            }
        }
    }

    // set the source
    d_source->set(srcids, shapes, strengths);
}

} // end namespace spn

//---------------------------------------------------------------------------//
//                 end of Problem_Builder.cc
//---------------------------------------------------------------------------//
