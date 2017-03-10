//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Work_Pool.cuh
 * \author Steven Hamilton
 * \date   Wed Mar 08 11:15:27 2017
 * \brief  Work_Pool class declaration.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Work_Pool_cuh
#define CudaUtils_cuda_utils_Work_Pool_cuh

#include <thrust/device_vector.h>

#include "CudaDBC.hh"
#include "Device_Memory_Manager.hh"
#include "Device_View_Field.hh"
#include "Utility_Functions.hh"

namespace cuda
{

//===========================================================================//
/*!
 * \class Work_Pool
 * \brief Class for pooling and distributing work among threads in a warp.
 */
/*!
 * \example cuda_utils/test/tstWork_Pool.cc
 *
 * Test of Work_Pool.
 */
//===========================================================================//

class Work_Pool
{
  public:
    //@{
    //! Typedefs
    typedef Device_View_Field<int>  Int_View;
    typedef Device_View_Field<bool> Bool_View;
    //@}

  private:
    // >>> DATA
    Int_View  d_indices;
    Bool_View d_active;

    int d_warp_num_work;
    int d_warp_begin;

  public:

    // Constructor
    Work_Pool(Int_View  indices,
              Bool_View active)
      : d_indices(indices)
      , d_active(active)
    {
        REQUIRE(d_indices.size() == d_active.size());
    }

    __device__ void initialize()
    {
        // Determine number of threads and number of warps
        int num_threads = blockDim.x * gridDim.x;
        int num_warps   = num_threads / warpSize;
        if (num_threads % warpSize)
            num_warps++;

        // Get id of warp
        int tid = cuda::utility::thread_id();
        int warp_id = tid / warpSize;

        // Determine number of indices and offset for this warp
        d_warp_num_work = d_indices.size() / num_warps;
        int extra_work  = d_indices.size() % num_warps;
        d_warp_begin = warp_id * d_warp_num_work;
        if (warp_id < extra_work)
        {
            d_warp_begin += warp_id;
            d_warp_num_work++;
        }
        else
        {
            d_warp_begin += extra_work;
        }
    }

    // Consolidate active indices
    __device__ void consolidate()
    {
        // Perform warp scan on current active indices
        int lane = threadIdx.x % warpSize;
        int work_sum = 0;
        for (int i = 0; i < work_per_thread(); ++i)
        {
            int id = lane + i * warpSize;

            // Get thread-local value for scan
            int active = 0;
            if (id < d_warp_num_work)
                active = d_active[id+d_warp_begin];

            // Compute exclusive scan of thread-local values
            int val, tmp_sum;
            utility::warpScan<int>(active,val,tmp_sum);

            // Update indices
            if (active)
            {
                // Indices should only move forward in the vector
                DEVICE_CHECK(work_sum+val <= id);
                
                // This looks like a potential race condition, but is not
                // because threads in a warp are always synchronized
                d_indices[work_sum+val+d_warp_begin] =
                    d_indices[id+d_warp_begin];
            }
            work_sum += tmp_sum;
        }

        // Now shuffle active into a consistent ordering
        for (int i = 0; i < work_per_thread(); ++i)
        {
            int id = lane + i * warpSize;
            if (id < d_warp_num_work)
            {
                if (id < work_sum)
                    d_active[id+d_warp_begin] = true;
                else
                    d_active[id+d_warp_begin] = false;
            }
        }
        d_warp_num_work = work_sum;
    }

    // Recompute active indices
    __device__ void set_inactive(int index)
    {
        DEVICE_REQUIRE(index < work_per_thread());

        int lane = threadIdx.x % warpSize;
        int work_id = lane + index * warpSize;
        d_active[work_id+d_warp_begin] = false;
    }

    // Amount of work per warp
    __device__ int work_per_warp() const
    {
        return d_warp_num_work;
    }

    // Amount of work per thread
    __device__ int work_per_thread() const
    {
        int work = d_warp_num_work / warpSize;
        if (d_warp_num_work % warpSize)
            work++;
        return work;
    }

    // Get global index for this thread given a work index
    __device__ int work_id(int index) const
    {
        DEVICE_REQUIRE(index < work_per_thread());

        int lane = threadIdx.x % warpSize;
        int work_id = lane + index * warpSize;
        int id = -1;
        if (work_id < d_warp_num_work && d_active[work_id+d_warp_begin])
            id = d_indices[work_id+d_warp_begin];
        return id;
    }
};

//===========================================================================//
/*!
 * \class Work_Pool_DMM
 * \brief Device_Memory_Manager for Work_Pool class.
 */
//===========================================================================//

class Work_Pool_DMM : public Device_Memory_Manager<Work_Pool>
{
  private:

      // >>>DATA

      thrust::device_vector<int>  d_indices;
      thrust::device_vector<bool> d_active;

  public:

      Work_Pool_DMM(const thrust::device_vector<int> &indices)
      {
          d_indices.assign(indices.begin(),indices.end());
          d_active.resize(d_indices.size());
      }

      Work_Pool device_instance()
      {
          thrust::fill(d_active.begin(),d_active.end(),true);
          Work_Pool pool(cuda::make_view(d_indices),
                         cuda::make_view(d_active));
          return pool;
      }
};

//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
// #include "Work_Pool.i.cuh"
//---------------------------------------------------------------------------//
#endif // CudaUtils_cuda_utils_Work_Pool_cuh

//---------------------------------------------------------------------------//
// end of CudaUtils/cuda_utils/Work_Pool.cuh
//---------------------------------------------------------------------------//
