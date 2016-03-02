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
    typedef cuda::arch::Device                  Arch_t;
    typedef cuda::Device_Vector<Arch_t,double>  Dev_Dbl;
    typedef cuda::Device_Vector<Arch_t,int>     Dev_Int;
    //@}

  private:
    // >>> DATA

    // Geometry and physics.
    Geometry_t *d_geometry;
    Physics_t  *d_physics;

    // Vector of tally result (single value, only batch statistics)
    double *d_tally;

    // List of cells we're tallying
    int *d_cells;

    // Volumes of tally cells
    double *d_volumes;

    // Number of cells being tallied
    int d_num_cells;

  public:

    // Constructor.
    Cell_Tally(SDP_Geometry geometry, SDP_Physics physics);

    // Disallow copy and assignment
    Cell_Tally(const Cell_Tally &phys)            = delete;
    Cell_Tally& operator=(const Cell_Tally &phys) = delete;

    // Destructor
    ~Cell_Tally()
    {
        cudaFree(d_tally);
        cudaFree(d_cells);
        cudaFree(d_volumes);
    }

    // Set cell list for tally
    void set_cells(const std::vector<int> &cells);

    // Get tally results.
    std::vector<double> results() const
    {
        std::vector<double> host_tally(d_num_cells);
        cudaMemcpy( &host_tally[0], d_tally, d_num_cells*sizeof(double),
                    cudaMemcpyDeviceToHost );
        return host_tally;
    }

    // Do post-processing on first and second moments
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
