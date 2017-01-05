//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_Array_Tester.cu
 * \author Tom Evans
 * \date   Wed Jan 04 23:32:17 2017
 * \brief  RTK_Array_Tester member and kernel definitions.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/fill.h>

#include "RTK_Array_Tester.hh"
#include "../RTK_Cell.cuh"
#include "../RTK_Array.cuh"

//---------------------------------------------------------------------------//
// TYPES
//---------------------------------------------------------------------------//

using Core_Manager = cuda_profugus::RTK_Core_Array_DMM;
using Core_Array   = Core_Manager::Core_Array;
using Vector       = Core_Array::Space_Vector;
using State        = Core_Array::Geo_State_t;

//---------------------------------------------------------------------------//
// SIMPLECORE
//---------------------------------------------------------------------------//

void SimpleCore::run_test()
{
    Core_Manager dmm(*core);
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/RTK_Array_Tester.cu
//---------------------------------------------------------------------------//
