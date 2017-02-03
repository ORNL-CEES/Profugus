//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_Geometry_Tester.cu
 * \author Tom Evans
 * \date   Fri Feb 03 09:50:55 2017
 * \brief  RTK_Geometry_Tester member and kernel definitions.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RTK_Geometry_Tester.hh"
#include "../RTK_Geometry.cuh"

//---------------------------------------------------------------------------//
// TYPES
//---------------------------------------------------------------------------//

using Core_Geometry         = cuda_profugus::Core;
using Core_Geometry_Manager = cuda_profugus::Core_DMM;

//---------------------------------------------------------------------------//
// CORE
//---------------------------------------------------------------------------//

void Core::heuristic()
{
    // Build the manager
    Core_Geometry_Manager dmm(*geometry);

    // Get the host object
    auto device_geo = dmm.device_instance();
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/RTK_Geometry_Tester.cu
//---------------------------------------------------------------------------//
