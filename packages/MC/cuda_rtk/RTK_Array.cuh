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

#include <memory>
#include <thrust/device_vector.h>

#include "Utils/utils/Definitions.hh"
#include "CudaUtils/cuda_utils/Device_Memory_Manager.hh"
#include "CudaUtils/cuda_utils/Device_View.hh"
#include "CudaUtils/cuda_utils/CudaDBC.hh"
#include "CudaUtils/cuda_utils/Device_View_Field.hh"
#include "CudaUtils/cuda_utils/Device_Vector_Lite.hh"
#include "CudaUtils/cuda_utils/Utility_Functions.hh"
#include "MC/geometry/RTK_Cell.hh"
#include "MC/geometry/RTK_Array.hh"
#include "RTK_Cell.cuh"
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
    using Dim_Vector   = cuda_utils::Device_Vector_Lite<int, 3>;
    //@}

  private:
    // Types
    using Object_t    = T;
    using Object_View = cuda::const_Device_View_Field<Object_t>;
    using View_Int    = cuda::const_Device_View_Field<int>;
    using View_Dbl    = cuda::const_Device_View_Field<double>;

  private:
    // >>> DATA

    // Array sizes
    Dim_Vector d_N;

    // Layout of objects in array.
    View_Int d_layout;

    // Array of core objects.
    Object_View d_objects;

    // Array dimensions.
    View_Dbl d_x, d_y, d_z;

    // Corner and length.
    Space_Vector d_corner;
    Space_Vector d_length;

    // Level
    int d_level;

    // Reflecting faces
    cuda_utils::Device_Vector_Lite<int, 6> d_reflect;

  public:
    // Constructor.
    RTK_Array(int level, Dim_Vector N, View_Int layout, Object_View objects,
              View_Dbl x, View_Dbl y, View_Dbl z, Space_Vector corner,
              Space_Vector length);

    // >>> DEVICE FUNCTIONS

    // >>> TRACKING FUNCTIONS

    // Initialize a state.
    __device__
    inline void initialize(const Space_Vector &r, Geo_State_t &state) const;

    // Track to next boundary.
    __device__
    inline void distance_to_boundary(const Space_Vector &r,
                                     const Space_Vector &omega,
                                     Geo_State_t &state) const;

    // Update a state at collision sites.
    __device__
    inline void update_state(Geo_State_t &state) const;

    // Cross a surface.
    __device__
    inline void cross_surface(const Space_Vector &r, Geo_State_t &state);

    // >>> ACCESSORS

    // Level of array.
    __device__
    int level() const { return d_level; }

    // Find the object a point is in.
    __device__
    inline int find_object(const Space_Vector &r, Geo_State_t &state) const;
    __device__
    inline int find_object_on_boundary(const Space_Vector &r, int face,
                                       int face_type, Geo_State_t &state) const;

    // Array indexing.
    __device__
    inline int index(int i, int j, int k) const;

    // Return the current material id.
    __device__
    inline int matid(const Geo_State_t &state) const;

  private:
    // >>> IMPLEMENTATION

    // Transform to object coordinate system.
    __device__
    inline Space_Vector transform(const Space_Vector &r,
                                  const Geo_State_t &state) const;

    // Get object.
    __device__
    inline const Object_t& object(const Geo_State_t &state) const;

    // Determine boundary crossings at each level in the array.
    __device__
    inline void determine_boundary_crossings(Geo_State_t &state) const;

    // Update the coordinates of each level in the array.
    __device__
    inline void update_coordinates(const Space_Vector &r,
                                   Geo_State_t &state) const;

    // Cross surface into next array element.
    __device__
    inline void calc_high_face(Geo_State_t &state, int face_type,
                               int exiting_face) const;
    __device__
    inline void calc_low_face(Geo_State_t &state, int face_type,
                              int exiting_face) const;
};

//---------------------------------------------------------------------------//
// SUPPORTED ARRAY TYPES
//---------------------------------------------------------------------------//

//! Lattice array (array of RTK_Cell objects).
using Lattice_Array = RTK_Array<RTK_Cell>;

//! Core array (array of Lattice_Array objects).
using Core_Array = RTK_Array<Lattice_Array>;

//===========================================================================//
/*!
 * \class RTK_Lattice_Array_DMM
 * \brief Device_Memory_Manager for a Lattice_Array.
 */
//===========================================================================//

class RTK_Lattice_Array_DMM
    : public cuda::Device_Memory_Manager<Lattice_Array>
{
  public:
    //! Host type.
    using Host_Lattice_Array = profugus::RTK_Array<profugus::RTK_Cell>;

  private:
    // Types.
    using Base         = cuda::Device_Memory_Manager<Lattice_Array>;
    using DV           = cuda::Device_View<RTK_Cell>;
    using SP_DV        = std::shared_ptr<DV>;
    using Space_Vector = Lattice_Array::Space_Vector;
    using Dim_Vector   = Lattice_Array::Dim_Vector;

  private:
    // >>> DEVICE MEMORY MANAGEMENT DATA

    // Objects stored in this array.
    SP_DV d_objects;

    // Array dimensions.
    thrust::device_vector<double> d_x, d_y, d_z;

    // Array layout.
    thrust::device_vector<int> d_layout;

  public:
    // Constructor.
    RTK_Lattice_Array_DMM(const Host_Lattice_Array &host_array);

    // >>> DERIVED INTERFACE

    // Create a device instance.
    Lattice_Array device_instance();

  private:
    // >>> IMPLEMENTATION

    // Array sizes
    Dim_Vector d_N;

    // Coordinates in real space of lower-left coordinate of array.
    Space_Vector d_corner;

    // Length in each dimension.
    Space_Vector d_length;
};

//===========================================================================//
/*!
 * \class RTK_Core_Array_DMM
 * \brief Device_Memory_Manager for a Core_Array.
 *
 * \note The construction of RTK arrays is not necessarily memory-optimized as
 * the same RTK cells could exist in multiple arrays, but since each array is
 * getting built individually these become \b copies on the device.  However,
 * this may not be a realistic issue as a depleted fuel core would not really
 * have the same pin in multiple arrays anyway.  Additionally, any methods
 * that check for already-built RTK_Cell objects would be potentially very
 * expensive (expensive hashing of the RTK_Cell data).
 */
//===========================================================================//
#if 0
class RTK_Core_Array_DMM
    : public cuda::Device_Memory_Manager<Core_Array>
{
  public:
    //! Host type.
    using Host_Lattice_Array = profugus::RTK_Array<profugus::RTK_Cell>;
    using Host_Core_Array    = profugus::RTK_Array<Host_Lattice_Array>;

  private:
    // Types.
    using Base  = cuda::Device_Memory_Manager<Core_Array>;
    using SP_DV = std::shared_ptr< cuda::Device_View<Lattice_Array> >;

  private:
    // >>> DEVICE MEMORY MANAGEMENT DATA

    // Objects stored in this array.
    SP_DV d_objects;

  public:
    // Constructor.
    RTK_Core_Array_DMM(const Host_Core_Array &host_array);

    // >>> DERIVED INTERFACE

    // Create a device instance.
    Core_Array device_instance();

  private:
    // >>> IMPLEMENTATION

};
#endif
} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// CUDA FUNCTIONS
//---------------------------------------------------------------------------//

#include "RTK_Array.i.cuh"

#endif // MC_cuda_rtk_RTK_Array_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Array.cuh
//---------------------------------------------------------------------------//
