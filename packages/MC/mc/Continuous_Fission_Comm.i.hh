//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_physics/Continuous_Fission_Comm.i.hh
 * \author Thomas M. Evans
 * \date   Fri Jul 29 17:33:14 2011
 * \brief  Continuous_Fission_Comm template member definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_physics_Continuous_Fission_Comm_i_hh
#define mc_physics_Continuous_Fission_Comm_i_hh

namespace shift
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Physics>
Continuous_Fission_Comm<Physics>::Continuous_Fission_Comm(
    SP_Boundary_Mesh bnd_mesh)
    : d_bnd_mesh(bnd_mesh)
    , d_size_fs(Physics::fission_site_bytes())
    , d_lox(d_bnd_mesh->low_block_edge(def::X))
    , d_hix(d_bnd_mesh->high_block_edge(def::X))
    , d_loy(d_bnd_mesh->low_block_edge(def::Y))
    , d_hiy(d_bnd_mesh->high_block_edge(def::Y))
    , d_loz(d_bnd_mesh->low_block_edge(def::Z))
    , d_hiz(d_bnd_mesh->high_block_edge(def::Z))
{
    using def::X; using def::Y; using def::Z;

    Require (d_bnd_mesh);
    Require (d_size_fs);

    // return if overlap is not on
    if (!d_bnd_mesh->overlap()) return;

    // make the adjacency matrix for connectivity of neighbor blocks
    d_M.resize(27);
    std::fill(d_M.begin(), d_M.end(), -1);

    // max, min, and offset for the neighbors surrounding this block
    Dims min, max;

    // get offset, min, and max by dimension
    for (int d = 0; d < 3; ++d)
    {
        // minimum by dimension
        if (d_bnd_mesh->logical(d) > 0)
            min[d] = 0;
        else
            min[d] = 1;

        // maximum by dimension
        if (d_bnd_mesh->logical(d) < d_bnd_mesh->num_blocks(d) - 1)
            max[d] = 3;
        else
            max[d] = 2;
    }

    // logical block coordinates for a given neighbor
    int bi = 0, bj = 0, bk = 0;

    // loop through neighbors and count and store them
    int ptr = 0, index = 0;
    for (int k = min[Z]; k < max[Z]; ++k)
    {
        // logical block coordinates in k for this neighbor
        bk = d_bnd_mesh->logical(Z) - 1 + k;
        Check (bk >= 0 && bk < d_bnd_mesh->num_blocks(Z));

        for (int j = min[Y]; j < max[Y]; ++j)
        {
            // logical block coordinates in j for this neighbor
            bj = d_bnd_mesh->logical(Y) - 1 + j;
            Check (bj >= 0 && bj < d_bnd_mesh->num_blocks(Y));

            for (int i = min[X]; i < max[X]; ++i)
            {
                // adjacency matrix index
                index = M_ijk(i, j, k);
                Check (index >= 0 && index < 27);

                // do not include this domain (1,1,1)
                if (index != 13)
                {
                    Check (index != M_ijk(1, 1, 1));

                    // logical block coordinates in i for this neighbor
                    bi = d_bnd_mesh->logical(X) - 1 + i;
                    Check (bi >= 0 && bi < d_bnd_mesh->num_blocks(X));

                    // assign the pointer into the block id array for this
                    // neighbor
                    d_M[index] = ptr;

                    // calculate the actual block index for this neighbor
                    d_block_id.push_back(d_bnd_mesh->convert(bi, bj, bk));

                    // advance the ptr
                    ptr++;
                }
            }
        }
    }

    Check (ptr == d_block_id.size());

    // >>> SET-CONSTANT communicator
    //     make communicator for all blocks within a set, ie. blocks 0,N in
    //     set 1 has a communicator, blocks 0,N in set 2 has a communicator,
    //     etc.
    nemesis::split(d_bnd_mesh->set(), 0, d_set_const_comm);

    // allocate receive handles
    d_handles.resize(ptr);

    // allocate send/receive buffers
    d_recv_size.resize(ptr, -1);
    d_send_buffers.resize(ptr);
    d_recv_buffers.resize(ptr);

    Remember (nemesis::set_internal_comm(d_set_const_comm));
    Ensure (nemesis::node()  == d_bnd_mesh->block());
    Ensure (nemesis::nodes() == d_bnd_mesh->num_blocks());
    Remember (nemesis::reset_internal_comm());
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Gather fission sites back to their parent domain.
 */
template<class Physics>
void Continuous_Fission_Comm<Physics>::gather(
    const Physics            &physics,
    Fission_Site_Container_t &fission_sites)
{
    using def::X; using def::Y; using def::Z;

    Require (d_bnd_mesh->overlap());
    Require (d_handles.size() == d_block_id.size());
    Require (d_handles.size() == d_send_buffers.size());
    Require (d_handles.size() == d_recv_buffers.size());
    Require (d_handles.size() == d_recv_size.size());

    // first set all communications to intra-set
    nemesis::set_internal_comm(d_set_const_comm);
    Check (nemesis::nodes() == d_bnd_mesh->num_blocks());
    Check (nemesis::node() == d_bnd_mesh->block());

    // number of send/receives
    int num_msgs = d_handles.size();

    // post receives for size of fission site domains
    for (int r = 0; r < num_msgs; ++r)
    {
        Check (d_send_buffers[r].empty());
        Check (d_recv_buffers[r].empty());

        nemesis::receive_async(d_handles[r], &d_recv_size[r], 1,
                               d_block_id[r], 201);
    }

    // make a new local domain fission site container
    Fission_Site_Container_t local_container;

    // loop through fission sites on this domain and fill buffers for neighbor
    // domains
    Space_Vector r;
    int          nbor_id;
    for (fs_iterator s = fission_sites.begin(); s != fission_sites.end(); ++s)
    {
        // get the position of the fission site from the physics object
        r = physics.fission_site(*s);

        // determine the neighbor it needs to go to
        nbor_id = d_M[M_ijk(search(r[X], d_lox, d_hix),
                            search(r[Y], d_loy, d_hiy),
                            search(r[Z], d_loz, d_hiz))];

        // if this is in an overlap region, then put it into the buffer for
        // the neighbor cell
        if (nbor_id > -1)
        {
            Check (nbor_id >= 0 && nbor_id < num_msgs);

            // add this site to the send buffer for the neighbor
            d_send_buffers[nbor_id].push_back(*s);
        }

        // otherwise add it to the new local container
        else
        {
            Check (M_ijk(search(r[X], d_lox, d_hix),
                         search(r[Y], d_loy, d_hiy),
                         search(r[Z], d_loz, d_hiz)) == 13);
            local_container.push_back(*s);
        }
    }

    // send sizes (number of fission sites) to neighbors
    nemesis::Request send_handle;
    int              size;
    for (int p = 0; p < num_msgs; ++p)
    {
        size = d_send_buffers[p].size();
        nemesis::send_async(send_handle, &size, 1, d_block_id[p], 201);
        send_handle.wait();
    }

    // allocate the receive buffers on each node and post receives
    for (int r = 0; r < num_msgs; ++r)
    {
        // get the sizes
        d_handles[r].wait();
        Check (d_recv_size[r] >= 0);
        Check (!d_handles[r].inuse());

        // post a new receive for buffers that have sites
        if (d_recv_size[r] > 0)
        {
            // allocate space in the receive buffer
            d_recv_buffers[r].resize(d_recv_size[r]);

            // make a void * to the container
            void *buffer = &d_recv_buffers[r][0];

            // post the receive
            nemesis::receive_async(
                d_handles[r], reinterpret_cast<char *>(buffer),
                d_recv_size[r] * d_size_fs, d_block_id[r], 202);
        }
    }

    // send buffers out
    for (int s = 0; s < num_msgs; ++s)
    {
        Check (!send_handle.inuse());
        size = d_send_buffers[s].size();

        // only send out buffers with data in them
        if (size > 0)
        {
            // make a void * to the container
            const void *buffer = &d_send_buffers[s][0];

            // send the data
            nemesis::send_async(
                send_handle, reinterpret_cast<const char *>(buffer),
                size * d_size_fs, d_block_id[s], 202);
            send_handle.wait();

            // clear the send-buffer
            d_send_buffers[s].clear();
        }

        Check (d_send_buffers[s].empty());
    }

    // wait on the receive buffers and add the fission sites to the
    // local_container
    for (int r = 0; r < num_msgs; ++r)
    {
        Check (d_recv_size[r] ? d_handles[r].inuse() : !d_handles[r].inuse());

        // wait on buffers that have data
        if (d_recv_size[r] > 0)
        {
            Check (d_handles[r].inuse());

            // wait on the buffer
            d_handles[r].wait();

            // add sites to the local container
            for (fs_iterator s = d_recv_buffers[r].begin();
                 s != d_recv_buffers[r].end(); ++s)
            {
                local_container.push_back(*s);
            }

            // empty the receive buffer
            d_recv_buffers[r].clear();
        }

        Check (d_recv_buffers[r].empty());
    }

    // set back to regular communicator
    nemesis::reset_internal_comm();
    Check (nemesis::nodes() == d_bnd_mesh->num_domains());
    Check (nemesis::node() == d_bnd_mesh->domain());

    // do an O(1) swap of the fission containers returning the new local
    // container in the input fission sites container
    fission_sites.swap(local_container);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return info about a neighbor.
 *
 * This function returns information about a neighbor using relative
 * coordinates, e.g. for the \c z=1 plane
 * \verbatim
       |---------|---------|---------|
       |         |         |         |
    2  | (0,2,1) | (1,2,1) | (2,2,1) |
       |         |         |         |
       |---------|---------|---------|
       |         |         |         |
  y 1  | (0,1,1) | (1,1,1) | (2,1,1) | z = 1
       |         |         |         |
       |---------|---------|---------|
       |         |         |         |
    0  | (0,0,1) | (1,0,1) | (2,0,1) |
       |         |         |         |
       |---------|---------|---------|

            0         1         2

                      x
 * \endverbatim
 * where the \c (1,1,1) block is the current block.
 *
 * \param i relative i-position of neighbor in range \c [0,2]
 * \param j relative j-position of neighbor in range \c [0,2]
 * \param k relative k-position of neighbor in range \c [0,2]
 *
 * \return a Neighbor_Block_Info object that gives the block and domain ids
 * for the neighbor along with its \c (i,j,k) coordinates in the boundary
 * mesh; if there is no neighbor in that direction (for blocks adjacent to a
 * boundary), then the Neighbor_Block_Info::exists == false.
 */
template<class Physics>
typename Continuous_Fission_Comm<Physics>::Neighbor_Block_Info
Continuous_Fission_Comm<Physics>::neighbor_info(int i,
                                                int j,
                                                int k) const
{
    using def::X; using def::Y; using def::Z;

    Require (d_bnd_mesh);
    Require (d_bnd_mesh->overlap());
    Require (d_M.size() == 27);

    // make a neighbor block info object
    Neighbor_Block_Info info;

    // canonical index of neighbor (in range [0, 26])
    int index = M_ijk(i, j, k);
    Check (index >= 0 && index < 27);

    // if the neighbor exists populate it
    if (d_M[index] > -1)
    {
        Check (d_M[index] < d_block_id.size());

        info.exists    = true;
        info.block_id  = d_block_id[d_M[index]];
        info.domain_id = d_bnd_mesh->domain(info.block_id);

        // logical coordinates
        info.ijk[X] = d_bnd_mesh->logical(X) - 1 + i;
        info.ijk[Y] = d_bnd_mesh->logical(Y) - 1 + j;
        info.ijk[Z] = d_bnd_mesh->logical(Z) - 1 + k;

        Ensure (d_bnd_mesh->convert(info.ijk[X], info.ijk[Y], info.ijk[Z]) ==
                info.block_id);
    }

    // otherwise, set exists to false
    else
    {
        info.exists    = false;
        info.block_id  = 0;
        info.domain_id = 0;
    }

    return info;
}

} // end namespace shift

#endif // mc_physics_Continuous_Fission_Comm_i_hh

//---------------------------------------------------------------------------//
//                   end of mc_physics/Continuous_Fission_Comm.i.hh
//---------------------------------------------------------------------------//
