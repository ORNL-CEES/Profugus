//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_Array_Tester.cu
 * \author Tom Evans
 * \date   Wed Jan 04 23:32:17 2017
 * \brief  RTK_Array_Tester member and kernel definitions.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/fill.h>

#include "RTK_Array_Tester.hh"
#include "../RTK_Cell.cuh"
#include "../RTK_Array.cuh"

//---------------------------------------------------------------------------//
// TYPES
//---------------------------------------------------------------------------//

using Lattice_Manager = cuda_profugus::Lattice_Array_DMM;
using Lattice_Array   = cuda_profugus::Lattice_Array;
using Core_Manager    = cuda_profugus::Core_Array_DMM;
using Core_Array      = cuda_profugus::Core_Array;
using Vector          = Lattice_Array::Space_Vector;
using State           = Lattice_Array::Geo_State_t;

//---------------------------------------------------------------------------//
// SIMPLELATTICE
//---------------------------------------------------------------------------//

__global__
void lattice_kernel(
    Lattice_Array  array,
    int           *ints)
{
    State state;
    Vector r;
    int m = 0;

    r = {1.261, 2.44, 12.1};
    array.initialize(r, state);

    ints[++m] = state.level_coord[0][0];
    ints[++m] = state.level_coord[0][1];
    ints[++m] = state.level_coord[0][2];
    ints[++m] = state.region;
    ints[++m] = array.matid(state);

    r = {1.259, 1.27, 1.1};
    array.initialize(r, state);

    ints[++m] = state.level_coord[0][0];
    ints[++m] = state.level_coord[0][1];
    ints[++m] = state.level_coord[0][2];
    ints[++m] = state.region;
    ints[++m] = array.matid(state);

    r = {3.560000,   2.239887,   1.300000};
    array.initialize(r, state);

    ints[++m] = state.level_coord[0][0];
    ints[++m] = state.level_coord[0][1];
    ints[++m] = state.level_coord[0][2];
    ints[++m] = state.region;
    ints[++m] = array.matid(state);

    r = {1.570000,   0.931993,   2.700000};
    array.initialize(r, state);

    ints[++m] = state.level_coord[0][0];
    ints[++m] = state.level_coord[0][1];
    ints[++m] = state.level_coord[0][2];
    ints[++m] = state.region;
    ints[++m] = array.matid(state);

    r = {1.300000,   2.044919,   3.800000};
    array.initialize(r, state);

    ints[++m] = state.level_coord[0][0];
    ints[++m] = state.level_coord[0][1];
    ints[++m] = state.level_coord[0][2];
    ints[++m] = state.region;
    ints[++m] = array.matid(state);
}

//---------------------------------------------------------------------------//

void SimpleLattice::run_test()
{
    // Make DMM
    Lattice_Manager dmm(*lattice);

    // Get the host object
    auto array = dmm.device_instance();

    thrust::device_vector<int> ints(25, -1);

    lattice_kernel<<<1,1>>>(array, ints.data().get());

    int m = 0;

    thrust::host_vector<int> rints(ints.begin(), ints.end());

    EXPECT_EQ(1,  rints[++m]);
    EXPECT_EQ(1,  rints[++m]);
    EXPECT_EQ(0,  rints[++m]);
    EXPECT_EQ(1,  rints[++m]);
    EXPECT_EQ(10, rints[++m]);

    EXPECT_EQ(0,  rints[++m]);
    EXPECT_EQ(1,  rints[++m]);
    EXPECT_EQ(0,  rints[++m]);
    EXPECT_EQ(1,  rints[++m]);
    EXPECT_EQ(10, rints[++m]);

    EXPECT_EQ(2,  rints[++m]);
    EXPECT_EQ(1,  rints[++m]);
    EXPECT_EQ(0,  rints[++m]);
    EXPECT_EQ(0,  rints[++m]);
    EXPECT_EQ(5,  rints[++m]);

    EXPECT_EQ(1,  rints[++m]);
    EXPECT_EQ(0,  rints[++m]);
    EXPECT_EQ(0,  rints[++m]);
    EXPECT_EQ(0,  rints[++m]);
    EXPECT_EQ(3,  rints[++m]);

    EXPECT_EQ(1,  rints[++m]);
    EXPECT_EQ(1,  rints[++m]);
    EXPECT_EQ(0,  rints[++m]);
    EXPECT_EQ(0,  rints[++m]);
    EXPECT_EQ(0,  rints[++m]);
}

//---------------------------------------------------------------------------//
// SIMPLECORE
//---------------------------------------------------------------------------//

__global__
void core_kernel(
    Core_Array  array,
    int        *ints)
{
    State state;
    Vector r;
    int m = 0;

    r = {1.2, 0.2, 4.5};
    array.initialize(r, state);

    ints[++m] = state.level_coord[1][0];
    ints[++m] = state.level_coord[1][1];
    ints[++m] = state.level_coord[1][2];
    ints[++m] = state.level_coord[0][0];
    ints[++m] = state.level_coord[0][1];
    ints[++m] = state.level_coord[0][2];
    ints[++m] = state.region;

    r = {8.34, 2.3, -3.1};
    array.initialize(r, state);

    ints[++m] = state.level_coord[1][0];
    ints[++m] = state.level_coord[1][1];
    ints[++m] = state.level_coord[1][2];
    ints[++m] = state.level_coord[0][0];
    ints[++m] = state.level_coord[0][1];
    ints[++m] = state.level_coord[0][2];
    ints[++m] = state.region;
}

//---------------------------------------------------------------------------//

void SimpleCore::run_test()
{
    // Make DMM
    Core_Manager dmm(*core);

    // Get the host object
    auto array = dmm.device_instance();

    thrust::device_vector<int> ints(25, -1);
    core_kernel<<<1,1>>>(array, ints.data().get());

    int m = 0;

    thrust::host_vector<int> rints(ints.begin(), ints.end());

    EXPECT_EQ(0, rints[++m]);
    EXPECT_EQ(0, rints[++m]);
    EXPECT_EQ(2, rints[++m]);
    EXPECT_EQ(0, rints[++m]);
    EXPECT_EQ(0, rints[++m]);
    EXPECT_EQ(0, rints[++m]);
    EXPECT_EQ(1, rints[++m]);

    EXPECT_EQ(2, rints[++m]);
    EXPECT_EQ(0, rints[++m]);
    EXPECT_EQ(0, rints[++m]);
    EXPECT_EQ(0, rints[++m]);
    EXPECT_EQ(1, rints[++m]);
    EXPECT_EQ(0, rints[++m]);
    EXPECT_EQ(0, rints[++m]);
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/RTK_Array_Tester.cu
//---------------------------------------------------------------------------//
