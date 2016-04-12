//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Tallier.cuh
 * \author Thomas M. Evans and Seth Johnson
 * \date   Mon May 12 12:15:30 2014
 * \brief  Tallier class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Tallier_cuh
#define cuda_mc_Tallier_cuh

#include "cuda_utils/Shared_Device_Ptr.hh"
#include "Particle.cuh"
#include "Physics.cuh"
#include "Cell_Tally.cuh"
#include "Keff_Tally.cuh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Tallier
 * \brief Do tally operations.
 */
/*!
 * \example mc/test/tstTallier.cc
 *
 * Test of Tallier.
 */
//===========================================================================//

template <class Geometry>
class Tallier
{
  public:
    //@{
    //! Typedefs.
    typedef Geometry                                Geometry_t;
    typedef Physics<Geometry_t>                     Physics_t;
    typedef Particle<Geometry_t>                    Particle_t;
    typedef Cell_Tally<Geometry_t>                  Cell_Tally_t;
    typedef Keff_Tally<Geometry_t>                  Keff_Tally_t;
    typedef cuda::Shared_Device_Ptr<Geometry_t>     SDP_Geometry;
    typedef cuda::Shared_Device_Ptr<Physics_t>      SDP_Physics;
    typedef cuda::Shared_Device_Ptr<Cell_Tally_t>   SDP_Cell_Tally;
    typedef cuda::Shared_Device_Ptr<Keff_Tally_t>   SDP_Keff_Tally;
    //@}

  private:
    // >>> DATA

    // Geometry and physics.
    Geometry_t *d_geometry;
    Physics_t  *d_physics;

    // Pointer to each tally type for now
    // Need to figure out how to manage this in an extensible way
    Cell_Tally_t *d_cell_tally;
    Keff_Tally_t *d_keff_tally;

    // Host-side objects
    SDP_Cell_Tally d_cell_tally_host;
    SDP_Keff_Tally d_keff_tally_host;

  public:

    // Constructor.
    Tallier();

    // Set the geometry and physics classes.
    void set(SDP_Geometry geometry, SDP_Physics physics);

    // Add cell tally.
    void add_cell_tally(SDP_Cell_Tally tally);

    // Add keff tally.
    void add_keff_tally(SDP_Keff_Tally tally);

    // Swap with another tallier
    void swap(Tallier<Geometry_t> &rhs);

    // >>> HOST-SIDE TALLY OPERATIONS

    // Tell the tallies to begin active kcode cycles
    void begin_active_cycles();

    // Tell the tallies to begin a new cycle in a kcode calculation
    void begin_cycle();

    // Tell the tallies to end a cycle in a kcode calculation
    void end_cycle(double num_particles);

    // Finalize tallies.
    void finalize(double num_particles);

    // Reset tallies.
    void reset();

    // >>> DEVICE-SIDE TALLY OPERATIONS

    // Process path-length tally events.
    __device__ void path_length(double step, const Particle_t &p);

    // Tally any source events.
    __device__ void source(const Particle_t &p);

    // Perform all end-history tally tasks.
    __device__ void end_history();

};

//---------------------------------------------------------------------------//
// NON-MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Swap two talliers.
 *
 * This is useful for deactivating tallying (like in a KCODE problem).
 *
 * This provides a std-like swap solution using Koenig namespace lookup.
 */
template <class Geometry>
inline void swap(Tallier<Geometry> &a,
                 Tallier<Geometry> &b)
{
    a.swap(b);
}

} // end namespace cuda_mc

#include "Tallier.i.cuh"

#endif // cuda_mc_Tallier_cuh

//---------------------------------------------------------------------------//
//                 end of Tallier.cuh
//---------------------------------------------------------------------------//
