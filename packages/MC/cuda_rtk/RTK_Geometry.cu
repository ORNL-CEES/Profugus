//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_Geometry.cu
 * \author Tom Evans
 * \date   Mon Jan 30 00:14:27 2017
 * \brief  RTK_Geometry member and kernel definitions.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RTK_Geometry.cuh"

#include "Utils/harness/DBC.hh"

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
// RTK_GEOMETRY
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
RTK_Geometry::RTK_Geometry(
    Array_t array)
    : d_array(array)
{
}

//---------------------------------------------------------------------------//
// RTK_GEOMETRY_DMM
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
RTK_Geometry_DMM::RTK_Geometry_DMM(
    const Host_Geometry &geometry)
    : d_array_manager(geometry.array())
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a device instance.
 */
RTK_Geometry RTK_Geometry_DMM::device_instance()
{
    return RTK_Geometry(d_array_manager.device_instance());
}

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Geometry.cu
//---------------------------------------------------------------------------//
