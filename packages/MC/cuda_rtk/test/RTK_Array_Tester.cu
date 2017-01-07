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

using Lattice_Manager = cuda_profugus::RTK_Lattice_Array_DMM;
using Lattice_Array   = cuda_profugus::Lattice_Array;
// using Core_Manager    = cuda_profugus::RTK_Core_Array_DMM;
// using Core_Array      = cuda_profugus::Core_Array;
using Vector          = Lattice_Array::Space_Vector;
using State           = Lattice_Array::Geo_State_t;

//---------------------------------------------------------------------------//
// SIMPLELATTICE
//---------------------------------------------------------------------------//

void SimpleLattice::run_test()
{
}

//---------------------------------------------------------------------------//
// SIMPLECORE
//---------------------------------------------------------------------------//

void SimpleCore::run_test()
{
    // Core_Manager dmm(*core);
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/RTK_Array_Tester.cu
//---------------------------------------------------------------------------//
