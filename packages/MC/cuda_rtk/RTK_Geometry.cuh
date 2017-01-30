//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_Geometry.cuh
 * \author Tom Evans
 * \date   Mon Jan 30 00:14:27 2017
 * \brief  RTK_Geometry CUDA-device class declarations.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_RTK_Geometry_cuh
#define MC_cuda_rtk_RTK_Geometry_cuh

#include "MC/geometry/RTK_Geometry.hh"
#include "RTK_Array.cuh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class RTK_Geometry
 * \brief CUDA implementation of profugus::RTK_Geometry.
 *
 * For cuda implementations the RTK_Geometry is explicity built using
 * Core_Array.  Thus, we \b always have a 2-level geometry, which is
 * sufficient for PWR core-level problems and BWR lattices.
 */
/*!
 * \example cuda_rtk/test/tstRTK_Geometry_cuda.cc
 *
 * Test of RTK_Geometry.
 */
//===========================================================================//

class RTK_Geometry
{
  public:
    //@{
    using Array_t = Core_Array;
    //@}

  private:
    // >>> DATA

    // Underlying core array.
    Array_t d_array;

  public:
    // Constructor.
    RTK_Geometry(Array_t array);

    // >>> DEVICE FUNCTIONS
};

//===========================================================================//
/*!
 * \class RTK_Geometry_DMM
 * \brief Device_Memory_Manager for an RTK_Geometry
 */
//===========================================================================//

class RTK_Geometry_DMM
    : public cuda::Device_Memory_Manager<RTK_Geometry>
{
  private:
    // Types.
    using Base          = cuda::Device_Memory_Manager<RTK_Geometry>;
    using Array_DMM     = Core_Array_DMM;
    using Host_Geometry = profugus::Core;

  private:
    // >>> DEVICE MEMORY MANAGEMENT DATA

    // Memory management of underlying array
    Array_DMM d_array_manager;

  public:
    // Constructor.
    RTK_Geometry_DMM(const Host_Geometry &geometry);

    // >>> DERIVED INTERFACE

    // Create a device instance.
    RTK_Geometry device_instance();
};

//---------------------------------------------------------------------------//
// SUPPORTED GEOMETRIES
//---------------------------------------------------------------------------//

//! Geometry built on Core_Array.
using Core = RTK_Geometry;

//! Geometry memory manager for Core.
using Core_DMM = RTK_Geometry_DMM;

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// CUDA FUNCTIONS
//---------------------------------------------------------------------------//

#include "RTK_Geometry.i.cuh"

#endif // MC_cuda_rtk_RTK_Geometry_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Geometry.cuh
//---------------------------------------------------------------------------//
