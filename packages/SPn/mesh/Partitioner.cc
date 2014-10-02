//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Partitioner.cc
 * \author Thomas M. Evans
 * \date   Wed Feb 12 09:54:53 2014
 * \brief  Partitioner member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "harness/Warnings.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"
#include "utils/Container_Functions.hh"
#include "Partitioner.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param pl
 */
Partitioner::Partitioner(RCP_ParameterList pl)
    : d_domain(profugus::node())
    , d_nodes(profugus::nodes())
{
    using def::I; using def::J; using def::K;

    REQUIRE(!pl.is_null());

    // initialize the parameter list
    init_pl(pl);

    // set dimension
    d_dimension = pl->get<int>("dimension");
    INSIST(d_dimension == 2 || d_dimension == 3, "Dimension must be 2 or 3");

    // initialize blocks and sets
    d_num_sets  = pl->get<int>("num_sets");
    if( pl->isType<int>("num_blocks_i") )
    {
        d_Nb[I] = pl->get<int>("num_blocks_i");
        INSIST( pl->isType<int>("num_blocks_j"),
                "Must specify num_blocks_j if specifying num_blocks_i." );
        d_Nb[J] = pl->get<int>("num_blocks_j");
        d_num_blocks = d_Nb[I] * d_Nb[J];
    }
    else
    {
        // Compute "most square" decomposition
        REQUIRE( d_nodes % d_num_sets == 0 );
        d_num_blocks = d_nodes / d_num_sets;
        int blocks_i, blocks_j;
        blocks_i = std::sqrt(static_cast<double>(d_num_blocks));
        blocks_j = d_num_blocks / blocks_i;
        while( blocks_i * blocks_j != d_num_blocks )
        {
            TEUCHOS_ASSERT(blocks_i < d_num_blocks);
            blocks_i++;
            blocks_j = d_num_blocks / blocks_i;
        }
        d_Nb[I] = blocks_i;
        d_Nb[J] = blocks_j;
    }
    d_k_blocks   = pl->get<int>("num_z_blocks");
    CHECK(d_num_blocks == d_Nb[I] * d_Nb[J]);
    CHECK(d_num_blocks * d_num_sets == d_nodes);

    // build global edges of mesh
    CHECK(d_edges[I].empty() && d_edges[J].empty() && d_edges[K].empty());

    // determine whether the mesh grid is set with uniform spacings,
    // by supplying the cell-edges explicitly for each direction,
    // or nothing (allow for derivative classes to set d_edges)

    int num_cells = 0;
    double delta  = 0.0;

    if (pl->isParameter("num_cells_i"))
    {
        num_cells = pl->get<int>("num_cells_i");
        delta     = pl->get<double>("delta_x");
        VALIDATE(num_cells > 0, "num_cells_i must be postitive");
        VALIDATE(delta > 0., "delta_x must be postitive");

        build_uniform_edges(num_cells, delta, d_edges[I]);
    }
    else
    {
        const Array_Dbl &a = pl->get<Array_Dbl>("x_edges");
        d_edges[I].insert(d_edges[I].end(), a.begin(), a.end());
        VALIDATE(profugus::is_sorted(d_edges[I].begin(), d_edges[I].end()),
                 "Mesh edges along X axis are not monotonically increasing.");
    }

    if (pl->isParameter("num_cells_j"))
    {
        num_cells = pl->get<int>("num_cells_j");
        delta     = pl->get<double>("delta_y");
        VALIDATE(num_cells > 0, "num_cells_j must be postitive");
        VALIDATE(delta > 0., "delta_y must be postitive");

        build_uniform_edges(num_cells, delta, d_edges[J]);
    }
    else
    {
        const Array_Dbl &a = pl->get<Array_Dbl>("y_edges");
        d_edges[J].insert(d_edges[J].end(), a.begin(), a.end());
        VALIDATE(profugus::is_sorted(d_edges[J].begin(), d_edges[J].end()),
                 "Mesh edges along Y axis are not monotonically increasing.");
    }

    // process K if 3-dimensional
    if (d_dimension == 3)
    {
        if (pl->isParameter("num_cells_k"))
        {
            num_cells = pl->get<int>("num_cells_k");
            delta     = pl->get<double>("delta_z");
            INSIST(num_cells > 0, "num_cells_k must be postitive");
            INSIST(delta > 0., "delta_z must be postitive");

            build_uniform_edges(num_cells, delta, d_edges[K]);
        }
        else
        {
            const Array_Dbl &a = pl->get<Array_Dbl>("z_edges");
            d_edges[K].insert(d_edges[K].end(), a.begin(), a.end());
            VALIDATE(profugus::is_sorted(d_edges[K].begin(), d_edges[K].end()),
                     "Mesh edges along Z axis are not monotonically "
                     "increasing.");
        }
    }

    ENSURE(d_Nb[I] > 0);
    ENSURE(d_Nb[J] > 0);
    ENSURE(d_dimension > 0);
    ENSURE(d_edges[I].size() > 1);
    ENSURE(d_edges[J].size() > 1);
    ENSURE(dimension() == 3 ? d_edges[K].size() > 1 : d_edges[K].empty());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build and partition the mesh and indexer on each processor.
 *
 * The build function:
 *  - initializes d_block, d_k_blocks
 *  - creates d_mesh, d_indexer, d_data
 */
void Partitioner::build()
{
    using def::I; using def::J; using def::K;

    REQUIRE(d_num_sets > 0);
    REQUIRE(d_num_blocks > 0);
    REQUIRE(d_num_blocks * d_num_sets == d_nodes);
    REQUIRE(d_edges[I].size() > 1);
    REQUIRE(d_edges[J].size() > 1);
    REQUIRE(d_dimension == 3 ? d_edges[K].size() > 1 : d_edges[K].empty());

    // >>> SET LOCAL PARTITION DATA
    // determine block and set ids.
    d_set   = d_domain / d_num_blocks;
    d_block = d_domain - d_set * d_num_blocks;
    CHECK(d_set < d_num_sets);
    CHECK(d_block < d_num_blocks);

    // determine the I/J indices of this block/plane
    Dim_Vector index;
    index[J] = d_block / d_Nb[I];
    index[I] = d_block - index[J] * d_Nb[I];
    CHECK(index[J] < d_Nb[J]);
    CHECK(index[I] < d_Nb[I]);
    CHECK(index[I] + index[J] * d_Nb[I] == d_block);

    // >>> SET GLOBAL PARTITION DATA
    // the coordinates in each direction on this processor
    IJK_Vec_Dbl local_edges;

    // the number of cells per block in each direction
    IJ_Vec_Int global_num;

    set_spatial_partition(index, local_edges, global_num);

    CHECK(!local_edges[I].empty());
    CHECK(!local_edges[J].empty());
    CHECK(!local_edges[K].empty() || d_dimension == 2);
    CHECK(global_num[I].size() == d_Nb[I]);
    CHECK(global_num[J].size() == d_Nb[J]);

    // Set the number of z blocks
    if (d_dimension == 3)
    {
        // update number of z blocks if mesh does not divide evenly
        bool change_blocks = false;
        const size_type num_cells_k = d_edges[K].size() - 1;
        while (num_cells_k % d_k_blocks != 0)
        {
            change_blocks = true;
            --d_k_blocks;
        }

        if (change_blocks && d_domain == 0)
        {
            ADD_WARNING("Number of Z blocks is being reduced to "
                    << d_k_blocks);
        }
        CHECK(d_k_blocks > 0);
    }

    // >>> BUILD THE MESH
    if (d_dimension == 3)
    {
        d_mesh = Teuchos::rcp(
            new Mesh_t(local_edges[I], local_edges[J], local_edges[K], index[I],
                       index[J], d_k_blocks));
        CHECK(!d_mesh.is_null());
        CHECK(d_mesh->num_cells() == (local_edges[I].size() - 1)
               * (local_edges[J].size() - 1) * (local_edges[K].size() - 1));
    }
    else
    {
        d_mesh = Teuchos::rcp(
            new Mesh_t(local_edges[I], local_edges[J], index[I], index[J]));
        CHECK(!d_mesh.is_null());
        CHECK(d_mesh->num_cells() == (local_edges[I].size() - 1)
               * (local_edges[J].size() - 1));
    }

    // >>> BUILD GLOBAL MESH DATA
    if (d_dimension == 3)
    {
        d_data = Teuchos::rcp(
            new Global_Mesh_Data(d_edges[I], d_edges[J], d_edges[K]));
    }
    else
    {
        d_data = Teuchos::rcp(
            new Global_Mesh_Data(d_edges[I], d_edges[J]));
    }

    // >>> BUILD INDEXER
    d_indexer = Teuchos::rcp(
        new LG_Indexer(global_num[I], global_num[J], d_num_sets));

    ENSURE(!d_mesh.is_null());
    ENSURE(!d_data.is_null());
    ENSURE(!d_indexer.is_null());
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize the parameter list with defaults.
 */
void Partitioner::init_pl(RCP_ParameterList pl)
{
    // make a parameter list of defaults
    ParameterList defaults;

    // number of blocks
    defaults.set("num_z_blocks", 1);

    // number of sets
    defaults.set("num_sets", 1);

    // dimension
    defaults.set("dimension", 3);

    // update the parameter list
    pl->setParametersNotAlreadySet(defaults);

    // determine dimension by examing the existing list
    bool has_k = (pl->isParameter("num_cells_k") ||
                  pl->isParameter("z_edges"));
    if (!has_k)
    {
        // set dimension to 2
        pl->set("dimension", 2);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build global uniform cell edges.
 */
void Partitioner::build_uniform_edges(int      num_cells,
                                      double   delta,
                                      Vec_Dbl& edges)
{
    REQUIRE(num_cells > 0);
    REQUIRE(delta > 0.);

    edges.resize(num_cells + 1);
    for (size_t i = 0; i < edges.size(); ++i)
        edges[i] = delta * i;

    ENSURE(edges.size() == num_cells + 1);
    ENSURE(profugus::soft_equiv(edges.front(), 0.));
    ENSURE(profugus::soft_equiv(edges.back(), delta * num_cells));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set spatial partitioning data
 *
 * The default partitioning scheme
 */
void Partitioner::set_spatial_partition(const Dim_Vector& index,
                                        IJK_Vec_Dbl&      local_edges,
                                        IJ_Vec_Int&       global_num) const
{
    using def::I; using def::J; using def::K;

    REQUIRE(index[I] < d_Nb[I]);
    REQUIRE(index[J] < d_Nb[J]);

    // global number of cells in each direction
    Dim_Vector global_num_cells;
    global_num_cells[I] = d_edges[I].size() - 1;
    global_num_cells[J] = d_edges[J].size() - 1;
    global_num_cells[K] = (d_dimension == 3 ? d_edges[K].size() - 1 : 1);

    // the number of cells on the I/J directions on this processor
    Dim_Vector local_num_cells;

    // loop through I/J directions and determine base cells per proc
    for (int dir = 0; dir < K; ++dir)
    {
        CHECK(d_Nb[dir] > 0);

        // estimate the base cells per pe in this direction
        local_num_cells[dir] = global_num_cells[dir] / d_Nb[dir];
        CHECK(local_num_cells[dir] > 0);

        // get the number of cells that need to be appended in this direction
        size_type append_cells = global_num_cells[dir] % d_Nb[dir];
        CHECK(append_cells < d_Nb[dir]);

        // if Ip/Jp is less than the number of append cells then add 1 cell
        // (row/column) to this block
        if (index[dir] < append_cells)
            local_num_cells[dir]++;
    }
    // size the z coordinates
    local_num_cells[K] = global_num_cells[K];

    // Set global number of cells per block
    set_global_num(local_num_cells, global_num);

    // Set local edge coordinates
    set_local_edges(index[I], global_num[I], d_edges[I], local_edges[I]);
    set_local_edges(index[J], global_num[J], d_edges[J], local_edges[J]);

    // K edges are global, so just copy them
    local_edges[K].assign(d_edges[K].begin(), d_edges[K].end());

    ENSURE(local_edges[I].size() == local_num_cells[I] + 1);
    ENSURE(local_edges[J].size() == local_num_cells[J] + 1);
    ENSURE(local_edges[K].size() == local_num_cells[K] + 1 ||
            dimension() == 2);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Communicate the number of cells in each block along I/J
 */
void Partitioner::set_global_num(const Dim_Vector& local_num_cells,
                                 IJ_Vec_Int&       global_num) const
{
    using def::I; using def::J; using def::K;

    REQUIRE(local_num_cells[I] > 0);
    REQUIRE(local_num_cells[J] > 0);

    // global number of cells on each domain
    IJ_Vec_Int global_cells;
    global_cells[I].assign(d_num_blocks, 0);
    global_cells[J].assign(d_num_blocks, 0);

    for (int dir = 0; dir < K; ++dir)
    {
        // assign this processors cells in I/J to the global list
        global_cells[dir][d_block] = local_num_cells[dir];

        // fill up the global cell list for this direction
        profugus::global_max(&(global_cells[dir][0]), d_num_blocks);
    }

    // loop through the first row of I (J = 0) and calculate the numbers from
    // the global cell processor list as we did above to calculate the spatial
    // offset for the block
    global_num[I].assign(global_cells[I].begin(),
                         global_cells[I].begin() + d_Nb[I]);

    // loop through first column of J (I = 0) and calculate the numbers from
    // the global cell processor list as we did above to calculate the spatial
    // offset for the block
    global_num[J].resize(d_Nb[J]);
    for (int j = 0, n = 0; j < d_Nb[J]; ++j)
    {
        // calculate index into global cell processor list
        n = j * d_Nb[I];
        CHECK(n >= 0 && n < d_nodes);

        global_num[J][j] = global_cells[J][n];
    }

    ENSURE(global_num[I].size() == d_Nb[I]);
    ENSURE(global_num[J].size() == d_Nb[J]);
    ENSURE(global_num[I].front() > 0);
    ENSURE(global_num[J].front() > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build non-uniform cell edges for a processor.
 *
 * \param index Index of this block along the direction
 * \param global_num_cells Global number of cells in the direction
 * \param global_edges Global edges in this direction
 * \param local_edges Resulting local edges
 */
void Partitioner::set_local_edges(const unsigned int index,
                                  const Vec_Int&     global_num_cells,
                                  const Vec_Dbl&     global_edges,
                                  Vec_Dbl&           local_edges) const
{
    using def::I; using def::J; using def::K;

    REQUIRE(!global_num_cells.empty());
    REQUIRE(!global_edges.empty());
    REQUIRE(index < global_num_cells.size());

    int start = std::accumulate(
        global_num_cells.begin(), global_num_cells.begin() + index, 0);

    // loop over directions and calculate the coordinates of this mesh block
    // partition
    local_edges.assign(
        global_edges.begin() + start,
        global_edges.begin() + start + global_num_cells[index] + 1);

    ENSURE(!local_edges.empty());
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Partitioner.cc
//---------------------------------------------------------------------------//
