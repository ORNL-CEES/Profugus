//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_Cell_Tester.cu
 * \author Tom Evans
 * \date   Tue Nov 29 17:08:48 2016
 * \brief  RTK_Cell_Tester member and kernel definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RTK_Cell_Tester.hh"
#include "../RTK_Cell.cuh"

//---------------------------------------------------------------------------//
// TYPES
//---------------------------------------------------------------------------//

using Device_Cell = cuda_profugus::RTK_Cell;
using DMM         = cuda_profugus::RTK_Cell_DMM;

//---------------------------------------------------------------------------//
// SINGLE-SHELL
//---------------------------------------------------------------------------//

void Single_Shell::run()
{
    DMM d1(*pins[0]), d2(*pins[1]);

    auto p1 = d1.device_instance();
    auto p2 = d2.device_instance();
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/RTK_Cell_Tester.cu
//---------------------------------------------------------------------------//
