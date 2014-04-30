//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_physics/Continuous_Fission_Comm.hh
 * \author Thomas M. Evans
 * \date   Fri Jul 29 17:33:14 2011
 * \brief  Continuous_Fission_Comm class definition.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_physics_Continuous_Fission_Comm_hh
#define mc_physics_Continuous_Fission_Comm_hh

#include <algorithm>
#include <vector>

#include "harness/DBC.hh"
#include "comm/global.hh"
#include "comm/Request.hh"
#include "utils/Definitions.hh"
#include "utils/SP.hh"
#include "utils/Vector_Lite.hh"
#include "mc/Boundary_Mesh.hh"

namespace shift
{

//===========================================================================//
/*!
 * \class Continuous_Fission_Comm
 * \brief Communicates fission sites in overlapping regions to the "owned"
 * domain for continuous fission site containers.
 *
 * The physics class must provide the Physics interface documented in
 * Physics.  The Fission_Site_Container must provide the following
 * types/methods:
 * - \c size()
 * - \c begin()
 * - \c end()
 * - \c operator[]
 * - \c push_back()
 * - \c swap()
 * - \c resize()
 * - \c clear()
 * - \c Fission_Site_Container::const_iterator
 * .
 * All random-access stl-type containers can provide this interface
 * (ie. vector, deque).  Naturally, this class is a specialization for
 * continuous fission site containers.  Thus, the total size of the fission
 * container in bytes should be recoverable using:
 * \code
    Fission_Site_Container_t container;
    // ...

    size_t num_bytes = sizeof(Fission_Site_t) * container.size();
    void *buffer     = &container[0];
    nemesis::send(reinterpret_cast<const char *>(buffer),
                  num_bytes, node);
   \endcode
 *
 * Both \c Fission_Site_t and \c Fission_Site_Container_t are defined in the
 * Physics class. Since the Fission_Site_Container_t is a continuous container
 * with \c operator[] defined, dereferenceing the first item in the container
 * points to the beginning of the memory.
 *
 * This class only has utility when overlapping domains are defined in a
 * domain-decomposed parallel decomposition.  If no overlapping domains are
 * defined can_gather() returns false and gather() will throw a failure
 * assertion if called.
 *
 * \sa shift::MG_Physics
 */
/*!
 * \example mc_physics/test/tstFission_Comm.cc
 *
 * Test of Continuous_Fission_Comm.
 */
//===========================================================================//

template<class Physics>
class Continuous_Fission_Comm
{
  public:
    //@{
    //! Typedefs.
    typedef typename Physics::Fission_Site           Fission_Site_t;
    typedef typename Physics::Fission_Site_Container Fission_Site_Container_t;
    typedef mc::Boundary_Mesh                        Boundary_Mesh_t;
    typedef nemesis::SP<Boundary_Mesh_t>             SP_Boundary_Mesh;
    typedef typename Boundary_Mesh_t::Vec_Int        Vec_Int;
    typedef nemesis::Vector_Lite<int, 3>             Dims;
    typedef typename Physics::Space_Vector           Space_Vector;
    //@}

    //! Iterators for fission site container.
    typedef typename Fission_Site_Container_t::const_iterator fs_iterator;

    //! Neighbor block information.
    struct Neighbor_Block_Info
    {
        //! Does the block neighbor actually exist?
        bool exists;

        //! If it exists, then the block id is.
        int  block_id;

        //! The domain id is equal to the block id if the number of sets = 1.
        int  domain_id;

        //! Logical indices of neighbor in boundary mesh.
        Dims ijk;
    };

  private:
    // >>> DATA

    // Boundary mesh.
    SP_Boundary_Mesh d_bnd_mesh;

  public:
    // Constructor.
    Continuous_Fission_Comm(SP_Boundary_Mesh bnd_mesh);

    //! Communicator is needed only if overlap exists.
    bool can_gather() const { return d_bnd_mesh->overlap(); }

    // Gather fission sites to their parent domain.
    void gather(const Physics &physics, Fission_Site_Container_t &fission_sites);

    // >>> ACCESSORS

    //! Number of neighbors.
    int num_neighbors() const { return d_block_id.size(); }

    //! Get neighbor information.
    Neighbor_Block_Info neighbor_info(int i, int j, int k) const;

    //! Size of a single fission site in bytes.
    int fission_site_size() const { return d_size_fs; }

  private:
    // >>> DATA

    // Set-constant communicator for block-to-block communication within a
    // set.
    nemesis::Communicator_t d_set_const_comm;

    // Neighbor block id.
    Vec_Int d_block_id;

    // Adjacency matrix.
    Vec_Int d_M;

    // Asynchronous receive handles.
    std::vector<nemesis::Request> d_handles;

    // Size of a fission site.
    int d_size_fs;

    // Send/receive fission buffers.
    std::vector<Fission_Site_Container_t> d_send_buffers, d_recv_buffers;

    // Send/receive sizes.
    Vec_Int d_recv_size;

    // Edges of block without overlapping regions.
    double d_lox, d_hix, d_loy, d_hiy, d_loz, d_hiz;

    // Point searches.
    int search(double point, double low, double hi) const
    {
        if (point < low) return 0;
        else if (point > hi) return 2;
        return 1;
    }

    // Convert (i,j,k) to 0-27 for adjacency matrix.
    int M_ijk(int i, int j, int k) const
    {
        Require (i >= 0 && i < 3);
        Require (j >= 0 && j < 3);
        Require (k >= 0 && k < 3);
        return i + 3 * (j + 3 * k);
    }
};

} // end namespace shift

//---------------------------------------------------------------------------//
// TEMPLATE MEMBER DEFINITIONS FOR AUTOMATIC INSTANTIATION
//---------------------------------------------------------------------------//

#include "Continuous_Fission_Comm.i.hh"

#endif // mc_physics_Continuous_Fission_Comm_hh

//---------------------------------------------------------------------------//
//              end of mc_physics/Continuous_Fission_Comm.hh
//---------------------------------------------------------------------------//
