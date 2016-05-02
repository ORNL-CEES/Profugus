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
    int *d_cells;

    // Vector of tally result (single value per entry, only batch statistics)
    double *d_tally;

    // Number of cells being tallied
    int d_num_cells;

    // >>> HOST-SIDE DATA

    // Shared pointers
    SDP_Geometry d_geometry_shared;

    // Thrust vectors
    thrust::device_vector<int>    d_cell_vec;
    thrust::device_vector<double> d_tally_vec;

    // Storage for host tally
    std::vector<int>    d_host_cells;
    std::vector<double> d_host_tally;
    std::vector<double> d_host_volumes;

  public:

    // Constructor.
    Cell_Tally(SDP_Geometry geometry, SDP_Physics physics);

    // Disallow copy and assignment
    Cell_Tally(const Cell_Tally &phys)            = delete;
    Cell_Tally& operator=(const Cell_Tally &phys) = delete;

    // Set cell list for tally
    void set_cells(const std::vector<int> &cells);

    // Get tally results.
    const std::vector<double> & results() const
    {
        REQUIRE( d_host_tally.size() == d_num_cells );
        return d_host_tally;
    }

    // Do post-processing
    void finalize(double num_particles);

    // Clear/re-initialize all tally values between solves
    void reset();

    // >>> ON-DEVICE FUNCTIONS

    // Accumulate tally for this history
    __device__ void accumulate(double step, const Particle_t &p);

    // Finalize for history
    // Null-op here because we're not computing a variance
    __device__ void end_history(){};

};

//---------------------------------------------------------------------------//

} // end namespace cuda_mc

#include "Cell_Tally.i.cuh"

#endif // MC_cuda_mc_Cell_Tally_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Cell_Tally.cuh
//---------------------------------------------------------------------------//
