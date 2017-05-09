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
#include "Particle_Vector.cuh"
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
    typedef Geometry                    Geometry_t;
    typedef Physics<Geometry_t>         Physics_t;
    typedef Particle_Vector<Geometry_t> Particle_Vector_t;
    typedef Cell_Tally<Geometry_t>      Cell_Tally_t;
    typedef Keff_Tally<Geometry_t>      Keff_Tally_t;
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

  public:

    // Constructor.
    Tallier(Geometry_t   *geometry,
            Physics_t    *physics,
            Cell_Tally_t *cell_tally,
            Keff_Tally_t *keff_tally)
        : d_geometry(geometry)
        , d_physics(physics)
        , d_cell_tally(cell_tally)
        , d_keff_tally(keff_tally)
    {
    }

    // >>> DEVICE-SIDE TALLY OPERATIONS

    // Process path-length tally events.
    __device__ void path_length(double step, int pid,
                                const Particle_Vector_t *particles);

    // Tally any source events.
    __device__ void source(int pid, const Particle_Vector_t *particles);

    // Perform all end-history tally tasks.
    __device__ void end_history();
};

//===========================================================================//
/*!
 * \class Tallier_DMM
 * \brief Device memory manager for Tallier class.
 */
//===========================================================================//

template <class Geometry>
class Tallier_DMM : public cuda::Device_Memory_Manager<Tallier<Geometry>>
{
  public:
    //@{
    //! Typedefs.
    typedef Geometry                                Geometry_t;
    typedef Physics<Geometry_t>                     Physics_t;
    typedef cuda::Shared_Device_Ptr<Geometry_t>     SDP_Geometry;
    typedef cuda::Shared_Device_Ptr<Physics_t>      SDP_Physics;
    typedef Cell_Tally<Geometry_t>                  Cell_Tally_t;
    typedef Keff_Tally<Geometry_t>                  Keff_Tally_t;
    typedef cuda::Shared_Device_Ptr<Cell_Tally_t>   SDP_Cell_Tally;
    typedef cuda::Shared_Device_Ptr<Keff_Tally_t>   SDP_Keff_Tally;
    typedef Cell_Tally_DMM<Geometry_t>              Cell_Tally_DMM_t;
    typedef Keff_Tally_DMM<Geometry_t>              Keff_Tally_DMM_t;
    typedef std::shared_ptr<Cell_Tally_DMM_t>       SP_Cell_Tally_DMM;
    typedef std::shared_ptr<Keff_Tally_DMM_t>       SP_Keff_Tally_DMM;
    //@}

  private:
    // >>> DATA

    // Pointer to each tally type for now
    // Need to figure out how to manage this in an extensible way
    SP_Cell_Tally_DMM   d_cell_tally_dmm;
    SP_Keff_Tally_DMM   d_keff_tally_dmm;
    SDP_Cell_Tally      d_cell_tally;
    SDP_Keff_Tally      d_keff_tally;

    // Host-side objects
    SDP_Geometry d_geometry;
    SDP_Physics  d_physics;

  public:

    // DMM interface
    Tallier<Geometry_t> device_instance()
    {
        if (d_cell_tally_dmm)
        {
            d_cell_tally = cuda::shared_device_ptr<Cell_Tally_t>(
                d_cell_tally_dmm->device_instance());
        }
        if (d_keff_tally_dmm)
        {
            d_keff_tally = cuda::shared_device_ptr<Keff_Tally_t>(
                d_keff_tally_dmm->device_instance());
        }
        Tallier<Geometry_t> tallier(d_geometry.get_device_ptr(),
                                    d_physics.get_device_ptr(),
                                    d_cell_tally.get_device_ptr(),
                                    d_keff_tally.get_device_ptr());
        return tallier;
    }

    // Set the geometry and physics classes.
    void set(SDP_Geometry geometry, SDP_Physics physics);

    // Add cell tally.
    void add_cell_tally(SP_Cell_Tally_DMM tally);

    // Add keff tally.
    void add_keff_tally(SP_Keff_Tally_DMM tally);

    // >>> HOST-SIDE TALLY OPERATIONS

    // Tell the tallies to begin active kcode cycles
    void begin_active_cycles();

    // Tell the tallies to begin a new cycle in a kcode calculation
    void begin_cycle(def::size_type num_particles);

    // Tell the tallies to end a cycle in a kcode calculation
    void end_cycle(double num_particles);

    // Tell the tallies to end a particle batch
    void end_batch(double num_particles);

    // Finalize tallies.
    void finalize(double num_particles);

    // Reset tallies.
    void reset();
};

} // end namespace cuda_mc

#include "Tallier.i.cuh"

#endif // cuda_mc_Tallier_cuh

//---------------------------------------------------------------------------//
//                 end of Tallier.cuh
//---------------------------------------------------------------------------//
