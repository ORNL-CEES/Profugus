//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/Linear_System_FV.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 14 19:58:19 2014
 * \brief  Linear_System_FV member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "comm/global.hh"
#include "comm/P_Stream.hh"
#include "utils/Constants.hh"

#include "Linear_System_FV.hh"

#include "Teuchos_DefaultComm.hpp"

namespace profugus
{
namespace tpetra
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR AND DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param dim SPN dimensions object
 * \param mat material database
 */
Linear_System_FV::Linear_System_FV(RCP_ParameterList db,
                                   RCP_Dimensions    dim,
                                   RCP_Mat_DB        mat,
                                   RCP_Mesh          mesh,
                                   RCP_Indexer       indexer,
                                   RCP_Global_Data   data,
                                   RCP_Timestep      dt)
    : Base(db, dim, mat, dt)
    , d_mesh(mesh)
    , d_indexer(indexer)
    , d_gather(mesh, b_mom_coeff, *indexer)
    , d_bnd_index(6)
    , d_Ng(mat->xs().num_groups())
    , d_Ne(dim->num_equations())
    , d_N(mesh->num_cells_dim(def::I),
          mesh->num_cells_dim(def::J),
          mesh->num_cells_dim(def::K))
    , d_Nc(d_N[def::I] * d_N[def::J] * d_N[def::K])
    , d_G(indexer->num_global(def::I),
          indexer->num_global(def::J),
          d_N[def::K])
    , d_Gc(d_G[def::I] * d_G[def::J] * d_G[def::K])
    , d_first(0)
    , d_last_I(d_G[def::I] - 1)
    , d_last_J(d_G[def::J] - 1)
    , d_last_K(d_G[def::K] - 1)
    , d_Nb_global(0)
    , d_Nb_local(0)
    , d_D_c(d_Ng, d_Ng)
    , d_C_c(d_Ng, d_Ng)
    , d_D(d_Ng, d_Ng)
    , d_C(d_Ng, d_Ng)
    , d_W(d_Ng, d_Ng)
    , d_widths(Vec_Dbl(data->num_cells(def::I)),
               Vec_Dbl(data->num_cells(def::J)),
               Vec_Dbl(data->num_cells(def::K)))
    , d_work(d_Ng)
    , d_ipiv(d_Ng)
{
    using def::I; using def::J; using def::K;

    REQUIRE(!mesh.is_null());
    REQUIRE(!indexer.is_null());
    REQUIRE(!data.is_null());
    REQUIRE(mesh->dimension() == 3);
    REQUIRE(data->which_dimension() == 3);
    REQUIRE(!b_dim.is_null());
    REQUIRE(!b_mat.is_null());
    REQUIRE(b_db->isParameter("boundary"));
    REQUIRE(d_Nc == d_mesh->num_cells());
    REQUIRE(d_Gc == data->num_cells());

    // make the communicator
    Teuchos::RCP<const Comm> comm = Teuchos::DefaultComm<int>::getComm();

    // only support single-set decompositions with SPN
    Insist(indexer->num_sets() == 1,
           "Only support 1-set decomposition in SPN.");

    // make the global cell widths
    for (int d = 0; d < def::END_IJK; ++d)
    {
        // get the global cell edges in this direction
        const Vec_Dbl &edges = data->edges(d);
        CHECK(edges.size() == d_widths[d].size() + 1);
        CHECK(edges.size() == data->num_cells(d) + 1);

        // loop through edges and generate the global cell width in each
        // "plane"
        for (int n = 0, N = d_widths[d].size(); n < N; ++n)
        {
            d_widths[d][n] = edges[n+1] - edges[n];
            CHECK(d_widths[d][n] > 0.0);
        }
    }

    // number of volume unknowns on this domain and global
    d_Nv_global = d_Gc * d_Ne * d_Ng;
    d_Nv_local  = d_Nc * d_Ne * d_Ng;

    // there are d_Ne * d_Ng unknowns per computational cell because energy is
    // not stored in a block matrix
    d_unknowns_per_cell = d_Ne * d_Ng;
    d_block_size        = 1;

    // calculate sizes of boundary data
    calc_bnd_sizes();

    // number of local and global rows in the matrix (number of cells X number
    // of equations + boundary unknowns)
    int N_local  = d_Nv_local  + d_Nb_local;
    int N_global = d_Nv_global + d_Nb_global;

    // make the block-map using the indexer; each block is an Ng x Ng matrix
    if (b_nodes > 1)
    {
        // local and global cell indices
        int global, local;

        // local-to-global cell map
        Vec_Int l2g(N_local, -1);

        // loop through cells on this mesh block and convert to global ids
        for (int k = 0; k < d_N[K]; ++k)
        {
            for (int j = 0; j < d_N[J]; ++j)
            {
                for (int i = 0; i < d_N[I]; ++i)
                {
                    // get the global cell index
                    global = d_indexer->l2g(i, j, k);
                    CHECK(global >= 0 && global < d_Gc);

                    // get the local cell index
                    local = d_mesh->convert(i, j, k);
                    CHECK(local >= 0 && local < d_Nc);

                    // loop over the equations (all equation moments are
                    // stored locally for a cell)
                    for (int n = 0; n < d_Ne; ++n)
                    {
                        // loop over all groups
                        for (int g = 0; g < d_Ng; ++g)
                        {
                            l2g[index(g, n, local)] = index(g, n, global);
                        }
                    }
                }
            }
        }

        // map boundary faces (these are appended after the volume cells) ->
        // these only exist when there are vacuum and/or source boundary
        // conditions
        map_bnd_l2g(d_bnd_index[0], d_N[J], d_N[K], l2g); // low  x face
        map_bnd_l2g(d_bnd_index[1], d_N[J], d_N[K], l2g); // high x face
        map_bnd_l2g(d_bnd_index[2], d_N[I], d_N[K], l2g); // low  y face
        map_bnd_l2g(d_bnd_index[3], d_N[I], d_N[K], l2g); // high y face
        map_bnd_l2g(d_bnd_index[4], d_N[I], d_N[J], l2g); // low  z face
        map_bnd_l2g(d_bnd_index[5], d_N[I], d_N[J], l2g); // high z face

        // make the map
        b_map = Teuchos::rcp(
            new Map_t(Teuchos::OrdinalTraits<int>::invalid(),
                           Teuchos::ArrayView<int>(l2g),0,comm) );
    }

    // otherwise we simply set the number of cells
    else
    {
        b_map = Teuchos::rcp(new Map_t(
                                    Teuchos::OrdinalTraits<int>::invalid(),
                                    N_global,0,comm) );
    }

    CHECK(b_map->getNodeNumElements() == N_local);
    CHECK(b_map->getGlobalNumElements() == N_global);

    // allocate space for the number of potential non-zero indices per row;
    // the matrix bandwidth (entries per row) is the (number of equations +
    // the number of spatially coupled cells) X the number of groups,
    // ie. (num-equations + 6) * Ng -> HOWEVER, we only ever insert at most 1
    // GXG block at a time, so this can be sized by G
    d_indices.resize(d_Ng);
    d_values.resize(d_Ng);

    // make the RHS vector
    b_rhs = Teuchos::rcp(new Vector_t(b_map));
    CHECK(b_rhs->getLocalLength() == N_local);
}

//---------------------------------------------------------------------------//
// INHERITED FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the matrix.
 */
void Linear_System_FV::build_Matrix()
{
    using def::I; using def::J; using def::K;

    REQUIRE(!d_mesh.is_null());
    REQUIRE(!d_indexer.is_null());

    // make the matrix
    d_matrix   = Teuchos::rcp(new CrsMatrix_t(b_map, d_Ng));
    b_operator = d_matrix;
    CHECK(!d_matrix->isStorageOptimized());

    // off-processor face fields of diffusion coefficients
    RCP_Face_Field Dx_low, Dx_high, Dy_low, Dy_high;

    // reference to indexer
    const LG_Indexer &index = *d_indexer;

    // global offsets in the (i,j) direction for this mesh block
    int i_off = index.offset(I);
    int j_off = index.offset(J);

    // global i,j and cell
    int g_i = 0, g_j = 0, global = 0;

    // >>> VOLUME EQUATIONS

    // outer loop over number of equations
    for (int eqn = 0; eqn < d_Ne; ++eqn)
    {
        // gather off-processor diffusion coefficients for this equation
        d_gather.gather(eqn);

        // get the face-fields holding the off-processor diffusion
        // coefficients (some may be null)
        Dx_low  = d_gather.low_side_D(I);
        Dy_low  = d_gather.low_side_D(J);
        Dx_high = d_gather.high_side_D(I);
        Dy_high = d_gather.high_side_D(J);

        for (int k = 0; k < d_N[K]; ++k)
        {
            for (int j = 0; j < d_N[J]; ++j)
            {
                for (int i = 0; i < d_N[I]; ++i)
                {
                    // get the global indices, we do not have to convert k
                    // because all of k lives on each processor for the KBA
                    // decomposition
                    g_i = i + i_off;
                    g_j = j + j_off;
                    CHECK(index.convert_to_global(i, j) ==
                           LG_Indexer::IJ_Set(g_i, g_j));
                    CHECK(index.l2g(i, j, k) == index.g2g(g_i, g_j, k));

                    // global cell index
                    global = index.l2g(i, j, k);

                    // build all of the G x G block matrices for this
                    // (equation, cell) block row; we will then add them
                    // element-by-element by inserting each row for each block
                    // matrix; this will incur some overhead, but everything
                    // will be summed and cleaned up at FillComplete()

                    // make the diffusion coefficient for this moment equation
                    // in this cell
                    b_mom_coeff->make_D(eqn, index.l2l(i, j, k), d_D_c);

                    // make A_nn matrix -> this adds A_nn to the
                    // diagonal-block (its placed in C_c); we need to do this
                    // before adding off-diagonal coupling terms
                    b_mom_coeff->make_A(eqn, eqn, index.l2l(i, j, k), d_C_c);

                    // FIRST: add spatially-coupled matrix elements
                    spatial_coupled_element(eqn, i, j, k, g_i, g_j,
                                            Dx_low, Dx_high, Dy_low, Dy_high);

                    // SECOND: insert the diagonal block
                    insert_block_matrix(eqn, global, 0, eqn, global, 0, d_C_c,
                                        *d_matrix);

                    // LAST: insert within-cell coupling with other moment
                    // equations
                    for (int m = 0; m < d_Ne; ++m)
                    {
                        if (m != eqn)
                        {
                            // make Anm
                            b_mom_coeff->make_A(eqn, m, index.l2l(i, j, k), d_W);

                            // add it to the matrix
                            insert_block_matrix(
                                eqn, global, 0, m, global, 0, d_W, *d_matrix);
                        }
                    }
                } // I
            } // J
        } // K
    } // eqn

    // >>> BOUNDARY EQUATIONS
    if (d_Nb_local)
    {
        // local and global face indices in each direction (-x,+x,-y,+y,-z,+z)
        int lf[6] = {0, d_N[I] - 1, 0, d_N[J] - 1, 0, d_N[K] - 1};
        int gf[6] = {d_first, d_last_I, d_first, d_last_J, d_first, d_last_K};

        // low/high x face
        for (int f = 0; f < 2; ++f)
        {
            // only process if this face has unknowns
            if (!d_bnd_index[f].is_null())
            {
                // loop over faces
                for (int k = 0; k < d_N[K]; ++k)
                {
                    for (int j = 0; j < d_N[J]; ++j)
                    {
                        // add boundary equation elements
                        add_boundary_equations(d_bnd_index[f]->l2g(j, k),
                                               d_indexer->l2g(lf[f], j, k),
                                               d_indexer->l2l(lf[f], j, k),
                                               d_widths[I][gf[f]]);
                    }
                }
            }
        }

        // low/high y face
        for (int f = 2; f < 4; ++f)
        {
            // only process if this face has unknowns
            if (!d_bnd_index[f].is_null())
            {
                // loop over faces
                for (int k = 0; k < d_N[K]; ++k)
                {
                    for (int i = 0; i < d_N[I]; ++i)
                    {
                        // add boundary equation elements
                        add_boundary_equations(d_bnd_index[f]->l2g(i, k),
                                               d_indexer->l2g(i, lf[f], k),
                                               d_indexer->l2l(i, lf[f], k),
                                               d_widths[J][gf[f]]);
                    }
                }
            }
        }

        // low/high z face
        for (int f = 4; f < 6; ++f)
        {
            // only process if this face has unknowns
            if (!d_bnd_index[f].is_null())
            {
                // loop over faces
                for (int j = 0; j < d_N[J]; ++j)
                {
                    for (int i = 0; i < d_N[I]; ++i)
                    {
                        // add boundary equation elements
                        add_boundary_equations(d_bnd_index[f]->l2g(i, j),
                                               d_indexer->l2g(i, j, lf[f]),
                                               d_indexer->l2l(i, j, lf[f]),
                                               d_widths[K][gf[f]]);
                    }
                }
            }
        }
    }

    // complete fill of matrix
    d_matrix->fillComplete();
    ENSURE(d_matrix->isStorageOptimized());

    // Epetra returns the global number of nonzeros as a 32 bit signed int
    //  which is prone to overflow (this is only used for output and doesn't
    //  represent any fundamental limitation).  We can get around this by
    //  global-summing the local number of nonzeros into a 64 bit int.
    UTILS_INT8 num_lhs_nonzeros = d_matrix->getNodeNumEntries();
    profugus::global_sum( num_lhs_nonzeros );

    profugus::pcout << ">>> Built SPN FV Element LHS Matrix with " <<
        num_lhs_nonzeros << " nonzero entries." << profugus::endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the right-hand-side fission matrix.
 */
void Linear_System_FV::build_fission_matrix()
{
    using def::I; using def::J; using def::K;

    REQUIRE(!d_mesh.is_null());
    REQUIRE(!b_mat.is_null());

    // as for the A matrix, here we build the matrix/graph at the same time
    // (since the coupling is much easier->no spatial coupling)
    d_fission = Teuchos::rcp(new CrsMatrix_t(b_map, d_Ng));
    b_fission = d_fission;

    // material id and (local/global) cell index
    int matid  = 0;
    int local  = 0;
    int global = 0;

    // outer loop over cells -> we can loop directly over cells because there
    // is no neighbor coupling in the fission matrix
    for (int k = 0; k < d_N[K]; ++k)
    {
        for (int j = 0; j < d_N[J]; ++j)
        {
            for (int i = 0; i < d_N[I]; ++i)
            {
                // get the cell indices
                global = d_indexer->l2g(i, j, k);
                local  = d_indexer->l2l(i, j, k);

                // inner loop over equations (elements in the row)
                for (int eqn = 0; eqn < d_Ne; ++eqn)
                {
                    // insert coupling with other moment equations
                    for (int m = 0; m < d_Ne; ++m)
                    {
                        // make Fnm
                        b_mom_coeff->make_F(eqn, m, local, d_W);

                        // add it to the matrix
                        insert_block_matrix(eqn, global, 0, m, global, 0,
                                            d_W, *d_fission);
                    }
                }
            }
        }
    }

    // finish matrix
    d_fission->fillComplete();

    ENSURE(d_fission->isLocallyIndexed());
    ENSURE(d_fission->isStorageOptimized());

    // Epetra returns the global number of nonzeros as a 32 bit signed int
    //  which is prone to overflow (this is only used for output and doesn't
    //  represent any fundamental limitation).  We can get around this by
    //  global-summing the local number of nonzeros into a 64 bit int.
    UTILS_INT8 num_rhs_nonzeros = d_fission->getNodeNumEntries();
    profugus::global_sum( num_rhs_nonzeros );

    profugus::pcout << ">>> Built SPN FV Element RHS Matrix with " <<
        num_rhs_nonzeros << " nonzero entries." << profugus::endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the RHS from an external (isotropic) source.
 *
 * Only isotropic external sources are allowed in this solver.
 *
 * \sa denovo::Source_DB
 */
void Linear_System_FV::build_RHS(const External_Source &q)
{
    using def::I; using def::J; using def::K;

    REQUIRE(q.num_groups() == d_Ng);
    REQUIRE(b_rhs->getLocalLength() == (d_Nv_local + d_Nb_local)
             * d_block_size);

    // reference to rhs vector
    Vector_t &rhs = *b_rhs;

    // set everything to zero
    rhs.putScalar(0.0);

    // >>> Add External Sources
    Insist(q.is_isotropic(),
            "Only isotropic sources are allowed in the SPN solver.");

    Teuchos::ArrayRCP<double> rhs_data = rhs.getDataNonConst(0);

    // loop through cells
    for (int cell = 0; cell < d_Nc; ++cell)
    {
        // loop over equations
        for (int n = 0; n < d_Ne; ++n)
        {
            // loop over groups
            for (int g = 0; g < d_Ng; ++g)
            {
                rhs_data[index(g, n, cell)] =
                    b_src_coefficients[n] * q.q_e(cell, g) *
                    profugus::constants::four_pi;
            }
        }
    }

    // if there are no boundary sources we are finished
    if (!d_Nb_local) return;

    // if the boundaries are vacuum we are also finished (sources are all
    // zero)
    if (b_db->get<std::string>("boundary") == "vacuum")
        return;

    // >>> Add Boundary Sources

    // local face index
    int face = d_Nc;

    // add boundary sources -> ORDER IS IMPORTANT HERE; we need to go in face
    // order, 0-1-2-3-4-5; since this is a private function I think it is ok
    // to allow this implicit convention that is formalized in the
    // FV_Bnd_indexer
    add_boundary_sources(0, face, "minus_x_phi", d_N[J]*d_N[K]); // low  x
    add_boundary_sources(1, face, "plus_x_phi",  d_N[J]*d_N[K]); // high x
    add_boundary_sources(2, face, "minus_y_phi", d_N[I]*d_N[K]); // low  y
    add_boundary_sources(3, face, "plus_y_phi",  d_N[I]*d_N[K]); // high y
    add_boundary_sources(4, face, "minus_z_phi", d_N[I]*d_N[J]); // low  z
    add_boundary_sources(5, face, "plus_z_phi",  d_N[I]*d_N[J]); // high z

    CHECK(face == d_Nc + d_Nb_local / d_unknowns_per_cell);
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate global and local boundary sizes.
 */
void Linear_System_FV::calc_bnd_sizes()
{
    using def::I; using def::J; using def::K;

    // get boundary condition info -> add space for boundary equations for any
    // non-reflecting faces

    // initialize all boundary faces to reflecting
    std::fill(d_bc_global, d_bc_global + 6, 0);
    std::fill(d_bc_local, d_bc_local + 6, 0);

    if (b_db->get<std::string>("boundary") == "reflect")
    {
        CHECK(b_db->isSublist("boundary_db"));

        // get the reflecting faces indicator
        const Array_Int &r = b_db->get<Teuchos::ParameterList>("boundary_db").
                             get<Array_Int>("reflect");
        CHECK(r.size() == 6);

        // set the number of boundary unknowns on each face
        d_bc_global[0] = (1 - r[0]) * d_G[J] * d_G[K];
        d_bc_global[1] = (1 - r[1]) * d_G[J] * d_G[K];
        d_bc_global[2] = (1 - r[2]) * d_G[I] * d_G[K];
        d_bc_global[3] = (1 - r[3]) * d_G[I] * d_G[K];
        d_bc_global[4] = (1 - r[4]) * d_G[I] * d_G[J];
        d_bc_global[5] = (1 - r[5]) * d_G[I] * d_G[J];
    }
    else
    {
        // we need to process edge fluxes on all faces
        d_bc_global[0] = d_G[J] * d_G[K];
        d_bc_global[1] = d_G[J] * d_G[K];
        d_bc_global[2] = d_G[I] * d_G[K];
        d_bc_global[3] = d_G[I] * d_G[K];
        d_bc_global[4] = d_G[I] * d_G[J];
        d_bc_global[5] = d_G[I] * d_G[J];
    }

    // calculate the global number of (block) boundary unknowns
    d_Nb_global = std::accumulate(d_bc_global, d_bc_global + 6, 0) *
                  d_unknowns_per_cell;

    // I/J partitions of mesh block
    int mesh_I = d_mesh->block(I);
    int mesh_J = d_mesh->block(J);

    // first and last block paritition indices
    int first  = 0;
    int last_I = d_indexer->num_blocks(I) - 1;
    int last_J = d_indexer->num_blocks(J) - 1;
    CHECK(mesh_I >= first && mesh_I <= last_I);
    CHECK(mesh_J >= first && mesh_J <= last_J);

    // this mesh block will have low I/J or high I/J boundary unknowns if it
    // is on the partition boundary and there are vaccum/source boundary
    // conditions on that boundary; every block has low/high K boundaries
    // because there is no paritioning in K

    // determine the local number of boundary faces

    // X-boundaries
    if (d_bc_global[0] > 0 && mesh_I == first)
    {
        d_bc_local[0] = d_N[J] * d_N[K];
    }
    if (d_bc_global[1] > 0 && mesh_I == last_I)
    {
        d_bc_local[1] = d_N[J] * d_N[K];
    }

    // Y-boundaries
    if (d_bc_global[2] > 0 && mesh_J == first)
    {
        d_bc_local[2] = d_N[I] * d_N[K];
    }
    if (d_bc_global[3] > 0 && mesh_J == last_J)
    {
        d_bc_local[3] = d_N[I] * d_N[K];
    }

    // Z-boundaries
    if (d_bc_global[4] > 0)
    {
        d_bc_local[4] = d_N[I] * d_N[J];
    }
    if (d_bc_global[5] > 0)
    {
        d_bc_local[5] = d_N[I] * d_N[J];
    }

    // calculate the local number of (block) boundary unknowns
    d_Nb_local = std::accumulate(d_bc_local, d_bc_local + 6, 0) *
                 d_unknowns_per_cell;

    // make local bnd face indexers; if there are local boundary cells on a
    // face it means (a) this is a boundary-adjacent block, and (b) there is a
    // vacuum or source boundary condition on the face -> NOTE: the offset for
    // K is always 0 since we are using the KBA decomposition
    if (d_bc_local[0])
    {
        CHECK(d_bc_global[0]);
        CHECK(mesh_I == first);
        d_bnd_index[0] = Teuchos::rcp(new FV_Bnd_Indexer(
            0, d_N[J], d_N[K], d_bc_local, d_G[J], d_G[K], d_bc_global,
            d_indexer->offset(J), 0));
    }
    if (d_bc_local[1])
    {
        CHECK(d_bc_global[1]);
        CHECK(mesh_I == last_I);
        d_bnd_index[1] = Teuchos::rcp(new FV_Bnd_Indexer(
            1, d_N[J], d_N[K], d_bc_local, d_G[J], d_G[K], d_bc_global,
            d_indexer->offset(J), 0));
    }
    if (d_bc_local[2])
    {
        CHECK(d_bc_global[2]);
        CHECK(mesh_J == first);
        d_bnd_index[2] = Teuchos::rcp(new FV_Bnd_Indexer(
            2, d_N[I], d_N[K], d_bc_local, d_G[I], d_G[K], d_bc_global,
            d_indexer->offset(I), 0));
    }
    if (d_bc_local[3])
    {
        CHECK(d_bc_global[3]);
        CHECK(mesh_J == last_J);
        d_bnd_index[3] = Teuchos::rcp(new FV_Bnd_Indexer(
            3, d_N[I], d_N[K], d_bc_local, d_G[I], d_G[K], d_bc_global,
            d_indexer->offset(I), 0));
    }
    if (d_bc_local[4])
    {
        CHECK(d_bc_global[4]);
        d_bnd_index[4] = Teuchos::rcp(new FV_Bnd_Indexer(
            4, d_N[I], d_N[J], d_bc_local, d_G[I], d_G[J], d_bc_global,
            d_indexer->offset(I), d_indexer->offset(J)));
    }
    if (d_bc_local[5])
    {
        CHECK(d_bc_global[5]);
        d_bnd_index[5] = Teuchos::rcp(new FV_Bnd_Indexer(
            5, d_N[I], d_N[J], d_bc_local, d_G[I], d_G[J], d_bc_global,
            d_indexer->offset(I), d_indexer->offset(J)));
    }

    ENSURE(d_Nb_local <= d_Nb_global);
    ENSURE(b_nodes == 1 ? d_Nb_local == d_Nb_global : true);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add boundary sources to the RHS vector
 */
void Linear_System_FV::add_boundary_sources(int         face_id,
                                            int        &face,
                                            const char *phi_d_str,
                                            int         num_face_cells)
{
    REQUIRE(face >= d_Nc);
    REQUIRE(b_db->isSublist("boundary_db"));

    // get a reference to the boundary sublist
    const Teuchos::ParameterList &bnd = b_db->sublist("boundary_db");

    // add sources on each face
    if (!d_bnd_index[face_id].is_null())
    {
        // check for sources on this face
        if (bnd.isParameter(phi_d_str))
        {
            CHECK(d_bc_global[face_id]); // better be sources on this face

            // reference to RHS vector
            Vector_t &rhs = *b_rhs;

            Teuchos::ArrayRCP<double> rhs_data = rhs.getDataNonConst(0);

            // get the boundary sources
            const Array_Dbl &phi_b = bnd.get<Array_Dbl>(phi_d_str);
            CHECK(phi_b.size() == d_Ng);

            // loop over faces
            for (int f = 0; f < num_face_cells; ++f, ++face)
            {
                // loop over equations
                for (int n = 0; n < d_Ne; ++n)
                {
                    // loop over groups
                    for (int g = 0; g < d_Ng; ++g)
                    {
                        rhs_data[index(g, n, face)] =
                            b_bnd_coefficients[n] * phi_b[g];
                    }
                }
            }
        }

        // advance the face index if there are no sources on this face
        // (ie. phi_b = 0.0 on this face), which will be the case for mixed
        // reflecting boundary conditions
        else
        {
            face += num_face_cells;
        }
    }

    ENSURE(face <= d_Nc + d_Nb_local / d_unknowns_per_cell);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map local-to-global boundary elements on a face.
 */
void Linear_System_FV::map_bnd_l2g(RCP_Bnd_Indexer  indexer,
                                   int              N_abscissa,
                                   int              N_ordinate,
                                   Vec_Int         &l2g)
{
    REQUIRE(l2g.size() == d_Nv_local + d_Nb_local);

    // only map elements on this boundary face if the indexer exists (if it
    // doesn't exist, then there are no boundary unknowns on this block)
    if (!indexer.is_null())
    {
        // loop over ordinate of face
        for (int o = 0; o < N_ordinate; ++o)
        {
            // loop over abscissa of face
            for (int a = 0; a < N_abscissa; ++a)
            {
                // loop over equations
                for (int n = 0; n < d_Ne; ++n)
                {
                    // loop over groups
                    for (int g = 0; g < d_Ng; ++g)
                    {
                        CHECK(d_Nv_local + index(g, n, indexer->local(a, o)) <
                               d_Nv_local + d_Nb_local);
                        CHECK(d_Nv_global + index(g, n, indexer->l2g(a, o)) <
                               d_Nv_global + d_Nb_global);

                        l2g[d_Nv_local + index(g, n, indexer->local(a, o))] =
                            d_Nv_global + index(g, n, indexer->l2g(a, o));
                    }
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Insert spatially coupled elements.
 */
void Linear_System_FV::spatial_coupled_element(int            n,
                                               int            i,
                                               int            j,
                                               int            k,
                                               int            g_i,
                                               int            g_j,
                                               RCP_Face_Field Dx_low,
                                               RCP_Face_Field Dx_high,
                                               RCP_Face_Field Dy_low,
                                               RCP_Face_Field Dy_high)
{
    using def::I; using def::J; using def::K;

    // global cell
    int global = d_indexer->l2g(i, j, k);
    CHECK(global == d_indexer->g2g(g_i, g_j, k));

    // spatially-coupled global cell index
    int neighbor = 0;

    if (g_i > d_first)
    {
        // neighbor cell
        neighbor = d_indexer->g2g(g_i - 1, g_j, k);

        // if we are on the inside of an internal block calculate the
        // diffusion coefficients locally
        if (i > 0)
        {
            // make the neighbor diffusion coefficient on this processor
            b_mom_coeff->make_D(n, d_indexer->l2l(i - 1, j, k), d_D);

            // add the spatial element to the matrix
            add_spatial_element(n, global, neighbor, d_widths[I][g_i - 1],
                                d_widths[I][g_i], d_widths[I][g_i], d_D_c, d_D);
        }
        // get the neighbor diffusion coefficient from the face-field if this
        // is on the low-internal-boundary side; we can only get here in a
        // multi-domain problem
        else
        {
            CHECK(!Dx_low.is_null());
            CHECK(b_nodes > 1);

            // this is a *View* into the field data, so we can't use a copy
            // matrix without potentially hosing data
            Serial_Matrix D = Dx_low->view(j, k);

            // add the spatial element to the matrix
            add_spatial_element(n, global, neighbor, d_widths[I][g_i - 1],
                                d_widths[I][g_i], d_widths[I][g_i], d_D_c, D);
        }
    }
    else if (!d_bnd_index[0].is_null())
    {
        // index of the edge-unknown
        neighbor = d_bnd_index[0]->l2g(j, k);

        // add the contribution from the boundary edge unknown
        build_bnd_element(n, global, neighbor, d_widths[I][d_first]);
    }

    if (g_i < d_last_I)
    {
        // neighbor cell
        neighbor = d_indexer->g2g(g_i + 1, g_j, k);

        // if we are on the inside of an internal block calculate the
        // diffusion coefficients locally
        if (i < d_N[I] - 1)
        {
            // make the neighbor diffusion coefficient on this processor
            b_mom_coeff->make_D(n, d_indexer->l2l(i + 1, j, k), d_D);

            // add the spatial element to the matrix
            add_spatial_element(n, global, neighbor, d_widths[I][g_i],
                                d_widths[I][g_i + 1], d_widths[I][g_i],
                                d_D, d_D_c);
        }
        // get the neighbor diffusion coefficient from the face-field if this
        // is on the high-internal-boundary side; we can only get here in a
        // multi-domain problem
        else
        {
            CHECK(!Dx_high.is_null());
            CHECK(b_nodes > 1);

            // this is a *View* into the field data, so we can't use a copy
            // matrix without potentially hosing data
            Serial_Matrix D = Dx_high->view(j, k);

            // add the spatial element to the matrix
            add_spatial_element(n, global, neighbor, d_widths[I][g_i],
                                d_widths[I][g_i + 1], d_widths[I][g_i],
                                D, d_D_c);
        }
    }
    else if (!d_bnd_index[1].is_null())
    {
        // index of the edge-unknown
        neighbor = d_bnd_index[1]->l2g(j, k);

        // add the contribution from the boundary edge unknown
        build_bnd_element(n, global, neighbor, d_widths[I][d_last_I]);
    }

    if (g_j > d_first)
    {
        // neighbor cell
        neighbor = d_indexer->g2g(g_i, g_j - 1, k);

        // if we are on the inside of an internal block calculate the
        // diffusion coefficients locally
        if (j > 0)
        {
            // make the neighbor diffusion coefficient on this processor
            b_mom_coeff->make_D(n, d_indexer->l2l(i, j - 1, k), d_D);

            // add the spatial element to the matrix
            add_spatial_element(n, global, neighbor, d_widths[J][g_j - 1],
                                d_widths[J][g_j], d_widths[J][g_j], d_D_c, d_D);
        }
        // get the neighbor diffusion coefficient from the face-field if this
        // is on the low-internal-boundary side; we can only get here in a
        // multi-domain problem
        else
        {
            CHECK(!Dy_low.is_null());
            CHECK(b_nodes > 1);

            // this is a *View* into the field data, so we can't use a copy
            // matrix without potentially hosing data
            Serial_Matrix D = Dy_low->view(i, k);

            // add the spatial element to the matrix
            add_spatial_element(n, global, neighbor, d_widths[J][g_j - 1],
                                d_widths[J][g_j], d_widths[J][g_j], d_D_c, D);
        }
    }
    else if (!d_bnd_index[2].is_null())
    {
        // index of the edge-unknown
        neighbor = d_bnd_index[2]->l2g(i, k);

        // add the contribution from the boundary edge unknown
        build_bnd_element(n, global, neighbor, d_widths[J][d_first]);
    }

    if (g_j < d_last_J)
    {
        // neighbor cell
        neighbor = d_indexer->g2g(g_i, g_j + 1, k);

        // if we are on the inside of an internal block calculate the
        // diffusion coefficients locally
        if (j < d_N[J] - 1)
        {
            // make the neighbor diffusion coefficient on this processor
            b_mom_coeff->make_D(n, d_indexer->l2l(i, j + 1, k), d_D);

            // add the spatial element to the matrix
            add_spatial_element(n, global, neighbor, d_widths[J][g_j],
                                d_widths[J][g_j + 1], d_widths[J][g_j],
                                d_D, d_D_c);
        }
        // get the neighbor diffusion coefficient from the face-field if this
        // is on the high-internal-boundary side; we can only get here in a
        // multi-domain problem
        else
        {
            CHECK(!Dy_high.is_null());
            CHECK(b_nodes > 1);

            // this is a *View* into the field data, so we can't use a copy
            // matrix without potentially hosing data
            Serial_Matrix D = Dy_high->view(i, k);

            // add the spatial element to the matrix
            add_spatial_element(n, global, neighbor, d_widths[J][g_j],
                                d_widths[J][g_j + 1], d_widths[J][g_j],
                                D, d_D_c);
        }
    }
    else if (!d_bnd_index[3].is_null())
    {
        // index of the edge-unknown
        neighbor = d_bnd_index[3]->l2g(i, k);

        // add the contribution from the boundary edge unknown
        build_bnd_element(n, global, neighbor, d_widths[J][d_last_J]);
    }

    if (k > d_first)
    {
        // neighbor cell
        neighbor = d_indexer->g2g(g_i, g_j, k - 1);

        // make the neighbor diffusion coefficient on this processor (K is
        // always local)
        b_mom_coeff->make_D(n, d_indexer->l2l(i, j, k - 1), d_D);
        add_spatial_element(n, global, neighbor, d_widths[K][k - 1],
                            d_widths[K][k], d_widths[K][k], d_D_c, d_D);
    }
    else if (!d_bnd_index[4].is_null())
    {
        // index of the edge-unknown
        neighbor = d_bnd_index[4]->l2g(i, j);

        // add the contribution from the boundary edge unknown
        build_bnd_element(n, global, neighbor, d_widths[K][d_first]);
    }

    if (k < d_last_K)
    {
        // neighbor cell
        neighbor = d_indexer->g2g(g_i, g_j, k + 1);

        // make the neighbor diffusion coefficient on this processor (K is
        // always local)
        b_mom_coeff->make_D(n, d_indexer->l2l(i, j, k + 1), d_D);
        add_spatial_element(n, global, neighbor, d_widths[K][k],
                            d_widths[K][k + 1], d_widths[K][k], d_D, d_D_c);
    }
    else if (!d_bnd_index[5].is_null())
    {
        // index of the edge-unknown
        neighbor = d_bnd_index[5]->l2g(i, j);

        // add the contribution from the boundary edge unknown
        build_bnd_element(n, global, neighbor, d_widths[K][d_last_K]);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief  Add spatial element to the matrix.
 *
 * \param eqn equation for this (GXG) block row
 * \param row_cell global cell for this (GXG) block row
 * \param col_cell global cell for this (GXG) block column
 * \param delta_l left-width
 * \param delta_r right-width
 * \param delta_c cell-width (will equal delta_l or delta_r)
 * \param Dl left-diffusion matrix
 * \param Dr right-diffusion matrix
 */
void Linear_System_FV::add_spatial_element(int                  eqn,
                                           int                  row_cell,
                                           int                  col_cell,
                                           double               delta_l,
                                           double               delta_r,
                                           double               delta_c,
                                           const Serial_Matrix &Dl,
                                           const Serial_Matrix &Dr)
{
    REQUIRE(Dl.numRows()    == d_Ng);
    REQUIRE(Dl.numCols()    == d_Ng);
    REQUIRE(Dr.numRows()    == d_Ng);
    REQUIRE(Dr.numCols()    == d_Ng);
    REQUIRE(d_C.numCols()   == d_Ng);
    REQUIRE(d_C.numRows()   == d_Ng);
    REQUIRE(d_C_c.numCols() == d_Ng);
    REQUIRE(d_C_c.numRows() == d_Ng);
    REQUIRE(delta_c > 0.0);

    // make C for this neighbor coupling -> note that Dl and Dr are references
    // to d_D and d_D_c depending on the spatial coupling direction

    // make the sum term (delta_l * Dl + delta_r * Dr)
    d_C.assign(Dl);
    d_C *= delta_l;

    d_W.assign(Dr);
    d_W *= delta_r;

    d_C += d_W;

    // invert

    // LU decomposition
    d_lapack.GETRF(d_Ng, d_Ng, d_C.values(), d_C.stride(), &d_ipiv[0], &d_info);
    CHECK(d_info == 0);

    // inverse
    d_lapack.GETRI(d_Ng, d_C.values(), d_C.stride(), &d_ipiv[0], &d_work[0],
                   d_Ng, &d_info);
    CHECK(d_info == 0);

    // multiply W = C*Dr
    d_W.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, d_C, Dr, 0.0);

    // multiply Dl * W = Dl * (C*Dr)
    d_C.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Dl, d_W, 0.0);

    // multiply by -2.0 / delta_c
    d_C *= -2.0 / delta_c;

    // add C to the spatially-coupled column
    insert_block_matrix(eqn, row_cell, 0,  eqn, col_cell, 0, d_C, *d_matrix);

    // add (C is negative) to the local (n,i,j,k) contribution to the matrix
    d_C_c -= d_C;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add boundary element to matrix.
 *
 * \param eqn
 * \param global_row
 * \param global_col
 * \param delta_c
 */
void Linear_System_FV::build_bnd_element(int    eqn,
                                         int    global_row,
                                         int    global_col,
                                         double delta_c)
{
    REQUIRE(d_C.numCols()   == d_Ng);
    REQUIRE(d_C.numRows()   == d_Ng);
    REQUIRE(d_C_c.numCols() == d_Ng);
    REQUIRE(d_C_c.numRows() == d_Ng);
    REQUIRE(delta_c > 0.0);

    // make C for this boundary coupling

    // make the sum term
    d_C.assign(d_D_c);
    d_C *= -2.0 / (delta_c * delta_c);

    // add C to the boundary-coupled column
    insert_block_matrix(eqn, global_row, 0, eqn, global_col, d_Nv_global, d_C,
                        *d_matrix);

    // add (C is negative) to the local (n,i,j,k) contribution to the matrix
    d_C_c -= d_C;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add boundary equations to the matrix.
 */
void Linear_System_FV::add_boundary_equations(int    face,
                                              int    global_cell,
                                              int    local_cell,
                                              double delta_c)
{
    // global row index for the current boundary equations
    int row = 0;

    // loop over equations
    for (int n = 0; n < d_Ne; ++n)
    {
        // make the diffusion coefficient for this moment equation in this
        // cell
        b_mom_coeff->make_D(n, local_cell, d_C);

        // calculate C
        d_C *= (-2.0 / delta_c);

        // make D_nn matrix and set it as the boundary term
        b_mom_coeff->make_B(n, n, d_C_c);

        // add C (C is negative) to the nn boundary term
        d_C_c -= d_C;

        // FIRST: add the volume cell term
        insert_block_matrix(n, face, d_Nv_global, n, global_cell, 0, d_C,
                            *d_matrix);

        // SECOND: add all the moment-coupling terms to the
        // boundary
        for (int m = 0; m < d_Ne; ++m)
        {
            if (m != n)
            {
                // make Bnm
                b_mom_coeff->make_B(n, m, d_W);

                // add it to the matrix
                insert_block_matrix(
                    n, face, d_Nv_global, m, face, d_Nv_global, d_W, *d_matrix);
            }
        }

        // LAST: add the within-moment matrix element
        insert_block_matrix(n, face, d_Nv_global, n, face, d_Nv_global, d_C_c,
                            *d_matrix);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Insert a block matrix (GXG) into the matrix.
 *
 * \param row_n equation for this (GXG) block row
 * \param row_cell global cell for this (GXG) block row
 * \param col_m equation for this (GXG) block column
 * \param col_cell global cell for this (GXG) block column
 * \param row_off used if these are boundary conditions that are appended to
 * the end of the matrix (should be d_Nv_global or 0)
 * \param col_off used if these are boundary conditions that are appended to
 * the end of the matrix (should be d_Nv_global or 0)
 */
void Linear_System_FV::insert_block_matrix(int                  row_n,
                                           int                  row_cell,
                                           int                  row_off,
                                           int                  col_m,
                                           int                  col_cell,
                                           int                  col_off,
                                           const Serial_Matrix &M,
                                           CrsMatrix_t         &matrix)
{
    REQUIRE(row_n < d_Ne);
    REQUIRE(col_m < d_Ne);
    REQUIRE(M.numCols() == d_Ng);
    REQUIRE(M.numRows() == d_Ng);

    // global row and column indices
    int row = 0, col = 0;

    // counter
    int ctr = 0;

    // loop over rows (in g)
    for(int g = 0; g < d_Ng; ++g)
    {
        // determine the row in global index space
        row = row_off + index(g, row_n, row_cell);
        CHECK(row >= 0 && row < d_Nv_global + d_Nb_global);

        // initialize the row counter
        ctr = 0;

        // loop over columns (in gp)
        for (int gp = 0; gp < d_Ng; ++gp)
        {
            CHECK(ctr < d_Ng);

            // determine the column in global index space
            col = col_off + index(gp, col_m, col_cell);
            CHECK(col >= 0 && col < d_Nv_global + d_Nb_global);

            // only add non-zero entries
            if (std::fabs(M(g, gp)) > 0.0)
            {
                d_indices[ctr] = col;
                d_values[ctr]  = M(g, gp);
                ++ctr;
            }
        }

        if( ctr >0 )
        {
            // insert the row into the matrix
            Teuchos::ArrayView<double> vals(&d_values[0],ctr);
            Teuchos::ArrayView<int> inds(&d_indices[0],ctr);
            matrix.insertGlobalValues(row, inds, vals );
        }
    }
}

} // end namespace profugus
} // end namespace tpetra

//---------------------------------------------------------------------------//
//                 end of Linear_System_FV.cc
//---------------------------------------------------------------------------//
