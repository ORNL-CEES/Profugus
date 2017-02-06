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
    using Object_t         = T;
    using Geo_State_t      = RTK_State;
    using Space_Vector     = Geo_State_t::Space_Vector;
    using Dim_Vector       = cuda_utils::Device_Vector_Lite<int, 3>;
    using Reflecting_Faces = cuda_utils::Device_Vector_Lite<int, 6>;
    //@}

  private:
    // Types
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
    Reflecting_Faces d_reflect;

  public:
    // Constructor.
    RTK_Array(int level, Dim_Vector N, View_Int layout, Object_View objects,
              View_Dbl x, View_Dbl y, View_Dbl z, Space_Vector corner,
              Space_Vector length, Reflecting_Faces reflect);

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

    // Find the object a point is in.
    __device__
    inline int find_object(const Space_Vector &r, Geo_State_t &state) const;
    __device__
    inline int find_object_on_boundary(const Space_Vector &r, int face,
                                       int face_type, Geo_State_t &state) const;

    // >>> ACCESSORS

    // Level of array.
    __host__ __device__
    int level() const { return d_level; }

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

//===========================================================================//
/*!
 * \class RTK_Array_DMM
 * \brief Device_Memory_Manager for an RTK_Array.
 */
//===========================================================================//

template<class Host_Array_T, class Device_Array_T>
class RTK_Array_DMM
    : public cuda::Device_Memory_Manager<Device_Array_T>
{
  private:
    // Types.
    using Base         = cuda::Device_Memory_Manager<Device_Array_T>;
    using Object_t     = typename Device_Array_T::Object_t;
    using DV           = cuda::Device_View<Object_t>;
    using SP_DV        = std::shared_ptr<DV>;
    using Space_Vector = typename Device_Array_T::Space_Vector;
    using Dim_Vector   = typename Device_Array_T::Dim_Vector;

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
    RTK_Array_DMM(const Host_Array_T &host_array);

    // >>> DERIVED INTERFACE

    // Create a device instance.
    Device_Array_T device_instance();

  private:
    // >>> IMPLEMENTATION

    // Array sizes
    Dim_Vector d_N;

    // Coordinates in real space of lower-left coordinate of array.
    Space_Vector d_corner;

    // Length in each dimension.
    Space_Vector d_length;

    // Level of the array.
    int d_level;

    // Reflecting faces.
    cuda_utils::Device_Vector_Lite<int, 6> d_reflect;
};

//---------------------------------------------------------------------------//
// SUPPORTED ARRAY TYPES
//---------------------------------------------------------------------------//

//! Lattice array (array of RTK_Cell objects).
using Lattice_Array = RTK_Array<RTK_Cell>;

//! Core array (array of Lattice_Array objects).
using Core_Array = RTK_Array<Lattice_Array>;

//! Host lattice array.
using Host_Lattice_Array = profugus::RTK_Array<profugus::RTK_Cell>;

//! Host core array.
using Host_Core_Array = profugus::RTK_Array<Host_Lattice_Array>;

//! Lattice array memory manager.
using Lattice_Array_DMM = RTK_Array_DMM<Host_Lattice_Array, Lattice_Array>;

//! Core array memory manager.
using Core_Array_DMM = RTK_Array_DMM<Host_Core_Array, Core_Array>;

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// CUDA FUNCTIONS
//---------------------------------------------------------------------------//

#include "RTK_Array.i.cuh"

#endif // MC_cuda_rtk_RTK_Array_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Array.cuh
//---------------------------------------------------------------------------//
