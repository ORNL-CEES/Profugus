//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/Cell_Tally.cuh
 * \author Thomas M. Evans
 * \date   Fri Aug 21 15:29:14 2015
 * \brief  Cell_Tally class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_mc_Cell_Tally_cuh
#define MC_cuda_mc_Cell_Tally_cuh

#include <vector>
#include <thrust/device_vector.h>

#include "cuda_utils/Shared_Device_Ptr.hh"
#include "cuda_utils/Device_Memory_Manager.hh"
#include "Physics.cuh"
#include "Particle.cuh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Cell_Tally
 * \brief Do pathlength cell tallies.
 */
/*!
 * \example mc/test/tstCell_Tally.cc
 *
 * Test of Cell_Tally.
 */
//===========================================================================//

template <class Geometry>
class Cell_Tally
{

  public:
    //@{
    //! Typedefs.
    typedef Geometry                            Geometry_t;
    typedef Physics<Geometry>                   Physics_t;
    typedef Particle<Geometry>                  Particle_t;
    typedef cuda::Shared_Device_Ptr<Geometry_t> SDP_Geometry;
    typedef cuda::Shared_Device_Ptr<Physics_t>  SDP_Physics;
    typedef thrust::device_vector<double>       Dev_Dbl;
    typedef thrust::device_vector<int>          Dev_Int;
    //@}

  private:

    // >>> DEVICE-SIDE DATA

    // Geometry and physics (device pointers)
    Geometry_t *d_geometry;
    Physics_t  *d_physics;

    // List of cells we're tallying
    const int *d_cells;

    // Vector of tally result (single value per entry, only batch statistics)
    double *d_tally;

    // Number of cells being tallied
    int d_num_cells;

  public:

    // Constructor.
    Cell_Tally(Geometry_t *geometry,
               Physics_t  *physics,
               const int  *cells,
               double     *tally,
               int         num_cells)
        : d_geometry(geometry)
        , d_physics(physics)
        , d_cells(cells)
        , d_tally(tally)
        , d_num_cells(num_cells)
    {
    }

    // >>> ON-DEVICE FUNCTIONS

    // Accumulate tally for this history
    __device__ void accumulate(double step, const Particle_t &p);

    // Finalize for history
    // Null-op here because we're not computing a variance
    __device__ void end_history(){};

};

//===========================================================================//
/*!
 * \class Cell_Tally_DMM
 * \brief Device memory manager for Cell_Tally
 */
//===========================================================================//

template <class Geometry>
class Cell_Tally_DMM : public cuda::Device_Memory_Manager<Cell_Tally<Geometry>>
{

  public:
    //@{
    //! Typedefs.
    typedef Geometry                            Geometry_t;
    typedef Physics<Geometry>                   Physics_t;
    typedef Particle<Geometry>                  Particle_t;
    typedef cuda::Shared_Device_Ptr<Geometry_t> SDP_Geometry;
    typedef cuda::Shared_Device_Ptr<Physics_t>  SDP_Physics;
    typedef thrust::device_vector<double>       Dev_Dbl;
    typedef thrust::device_vector<int>          Dev_Int;
    //@}

  private:


    // >>> HOST-SIDE DATA

    // Shared pointers
    SDP_Geometry d_geometry;
    SDP_Physics  d_physics;

    // Thrust vectors
    thrust::device_vector<int>    d_cells;
    thrust::device_vector<double> d_tally;

    // Storage for host tally
    std::vector<int>    d_host_cells;
    std::vector<double> d_host_tally;
    std::vector<double> d_host_std_dev;
    std::vector<double> d_host_volumes;

    // Data for batch statistics
    std::vector<double> d_tally_sum;
    std::vector<double> d_tally_sum_sq;
    std::vector<double> d_batch_np;

  public:

    // Constructor.
    Cell_Tally_DMM(SDP_Geometry geometry, SDP_Physics physics);

    // DMM interface
    Cell_Tally<Geometry> device_instance()
    {
        REQUIRE(d_cells.size()      == d_tally.size());
        REQUIRE(d_host_cells.size() == d_cells.size());

        Cell_Tally<Geometry> tally(d_geometry.get_device_ptr(),
                                   d_physics.get_device_ptr(),
                                   d_cells.data().get(),
                                   d_tally.data().get(),
                                   d_cells.size());
        return tally;
    }

    // Set cell list for tally
    void set_cells(const std::vector<int>    &cells,
                   const std::vector<double> &volumes);

    // Get tally results.
    const std::vector<double> & results() const
    {
        REQUIRE(d_host_tally.size() == d_cells.size());
        return d_host_tally;
    }

    // Get tally results.
    const std::vector<double> & std_dev() const
    {
        REQUIRE(d_host_std_dev.size() == d_cells.size());
        return d_host_std_dev;
    }

    // Finalize batch
    void end_batch(double num_particles);

    // Do post-processing
    void finalize(double num_particles);

    // Clear/re-initialize all tally values between solves
    void reset();
};

//---------------------------------------------------------------------------//

} // end namespace cuda_mc

#include "Cell_Tally.i.cuh"

#endif // MC_cuda_mc_Cell_Tally_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Cell_Tally.cuh
//---------------------------------------------------------------------------//
