//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_Array.cuh
 * \author Tom Evans
 * \date   Wed Jan 04 15:43:43 2017
 * \brief  RTK_Array CUDA-device class declarations.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_RTK_Array_cuh
#define MC_cuda_rtk_RTK_Array_cuh

#include <thrust/device_vector.h>

#include "CudaUtils/cuda_utils/Device_Memory_Manager.hh"
#include "CudaUtils/cuda_utils/Device_View.hh"
#include "CudaUtils/cuda_utils/CudaDBC.hh"
#include "CudaUtils/cuda_utils/Device_View_Field.hh"
#include "CudaUtils/cuda_utils/Device_Vector_Lite.hh"
#include "MC/geometry/RTK_Array.hh"
#include "RTK_State.cuh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class RTK_Array
 * \brief CUDA implementation of profugus::RTK_Array.
 */
/*!
 * \example cuda_rtk/test/tstRTK_Array_cuda.cc
 *
 * Test of RTK_Array.
 */
//===========================================================================//

template<class T>
class RTK_Array
{
  public:
    //@{
    //! Types.
    using Geo_State_t  = RTK_State;
    using Space_Vector = Geo_State_t::Space_Vector;
    //@}

  private:
    // Types
    using Object_View = cuda::const_Device_View_Field<T>;
    using View_Int    = cuda::const_Device_View_Field<int>;
    using View_Dbl    = cuda::const_Device_View_Field<double>;

  private:
    // >>> DATA

    // Array sizes
    cuda_utils::Device_Vector_Lite<int, 3> d_N;

    // Layout of objects in array.
    View_Int d_layout;

    // Array of core objects.
    Object_View d_objects;

    // Array dimensions.
    View_Dbl d_x, d_y, d_z;

  public:
    // Constructor.
    RTK_Array();

    // >>> DEVICE FUNCTIONS

    // Find the object a point is in.
    __device__
    int find_object(const Space_Vector &r, Geo_State_t &state) const;
};

//===========================================================================//
/*!
 * \class RTK_Array_DMM
 * \brief Device_Memory_Manager for the device-based RTK_Array class.
 *
 * \note The construction of RTK arrays is not necessarily memory-optimized as
 * the same RTK cells could exist in multiple arrays, but since each array is
 * getting built individually these become \b copies on the device.  However,
 * this may not be a realistic issue as a depleted fuel core would not really
 * have the same pin in multiple arrays anyway.
 */
//===========================================================================//

template<class T>
class RTK_Array_DMM
    : public cuda::Device_Memory_Manager< RTK_Array<T> >
{
    using Base = cuda::Device_Memory_Manager< RTK_Array<T> >;

  public:
    //! Host type.
    using Host_RTK_Array = profugus::RTK_Array<T>;

  private:
    // Types.
    using DV = cuda::Device_View<T>;

  private:
    // >>> DEVICE MEMORY MANAGEMENT DATA

    // Objects stored in this array.
    DV d_objects;

  public:
    // Constructor.
    RTK_Array_DMM(const Host_RTK_Array &host_array);

    // >>> DERIVED INTERFACE

    // Create a device instance.
    RTK_Array<T> device_instance();

  private:
    // >>> IMPLEMENTATION

};

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// CUDA FUNCTIONS
//---------------------------------------------------------------------------//

#include "RTK_Array.i.cuh"

#endif // MC_cuda_rtk_RTK_Array_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Array.cuh
//---------------------------------------------------------------------------//
