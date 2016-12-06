//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_Cell_Tester.cu
 * \author Tom Evans
 * \date   Tue Nov 29 17:08:48 2016
 * \brief  RTK_Cell_Tester member and kernel definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "RTK_Cell_Tester.hh"
#include "../RTK_Cell.cuh"

//---------------------------------------------------------------------------//
// TYPES
//---------------------------------------------------------------------------//

using Device_Cell = cuda_profugus::RTK_Cell;
using DMM         = cuda_profugus::RTK_Cell_DMM;
using Vector      = Device_Cell::Space_Vector;
using State       = Device_Cell::Geo_State_t;

using def::X; using def::Y; using def::Z;

//---------------------------------------------------------------------------//
// SINGLE-SHELL
//---------------------------------------------------------------------------//

__global__
void single_shell_kernel1(
    Device_Cell  pin,
    int         *ints,
    double      *dbls,
    Vector      *svs)
{
    // Test basic pin parameters
    {
        auto r  = pin.radii();
        dbls[0] = r[0];
        dbls[1] = pin.pitch(X);
        dbls[2] = pin.pitch(Y);

        ints[0] = pin.num_regions();
        ints[1] = pin.num_shells();
        ints[2] = pin.region(0.0, 0.5401);
        ints[3] = pin.region(0.0, 0.5399);
        ints[4] = pin.matid(0);
        ints[5] = pin.matid(1);

        pin.get_extents(svs[0], svs[1]);
    }

    // Make the state and initialize
    State state;
    {
        pin.initialize({0.0, 0.53, 0.0}, state);
        ints[6] = state.region;
        ints[7] = state.segment;
        pin.initialize({0.0, 0.54, 0.0}, state);
        ints[8]  = state.region;
        ints[9]  = state.segment;
        ints[10] = pin.cell(state.region, state.segment);
        pin.initialize({0.0, 0.55, 0.0}, state);
        ints[11] = state.region;
        ints[12] = state.segment;
        ints[13] = pin.cell(state.region, state.segment);
    }

    // Distance-to-boundary box-boundary tests
    Vector r, omega;
    {
        r     = {0.0, 0.59, 0.0};
        omega = {1.0, 0.0, 0.0};
        pin.initialize(r, state);
        pin.distance_to_boundary(r, omega, state);

        dbls[3]  = state.dist_to_next_region;
        ints[14] = state.exiting_face;
        ints[15] = state.region;
    }
    {
        r     = {0.0, 0.59, 0.0};
        omega = {-1.0, 0.0, 0.0};
        pin.initialize(r, state);
        pin.distance_to_boundary(r, omega, state);

        dbls[4]  = state.dist_to_next_region;
        ints[16] = state.exiting_face;
        ints[17] = state.region;
    }
}

//---------------------------------------------------------------------------//

__global__
void single_shell_kernel2(
    Device_Cell pin)
{
}

//---------------------------------------------------------------------------//

void Single_Shell::run()
{
    DMM d1(*pins[0]), d2(*pins[1]);

    auto p1 = d1.device_instance();
    auto p2 = d2.device_instance();

    thrust::device_vector<int>    ints(20, -1);
    thrust::device_vector<double> dbls(20, -1);
    thrust::device_vector<Vector> svs(20);

    {
        single_shell_kernel1<<<1,1>>>(
            p1, ints.data().get(), dbls.data().get(), svs.data().get());

        thrust::host_vector<int>    rints(ints.begin(), ints.end());
        thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());
        thrust::host_vector<Vector> rsvs(svs.begin(), svs.end());

        EXPECT_EQ(0.54, dbls[0]);
        EXPECT_EQ(1.26, dbls[1]);
        EXPECT_EQ(1.26, dbls[2]);

        EXPECT_EQ(2,  rints[0]);
        EXPECT_EQ(1,  rints[1]);
        EXPECT_EQ(1,  rints[2]);
        EXPECT_EQ(0,  rints[3]);
        EXPECT_EQ(1,  rints[4]);
        EXPECT_EQ(10, rints[5]);

        EXPECT_SOFTEQ(rsvs[0][X], -1.26/2, 1.e-12);
        EXPECT_SOFTEQ(rsvs[0][Y], -1.26/2, 1.e-12);
        EXPECT_SOFTEQ(rsvs[0][Z],  0.    , 1.e-12);
        EXPECT_SOFTEQ(rsvs[1][X],  1.26/2, 1.e-12);
        EXPECT_SOFTEQ(rsvs[1][Y],  1.26/2, 1.e-12);
        EXPECT_SOFTEQ(rsvs[1][Z], 14.28  , 1.e-12);

        EXPECT_EQ(0, rints[6]);
        EXPECT_EQ(0, rints[7]);
        EXPECT_EQ(0, rints[8]);
        EXPECT_EQ(0, rints[9]);
        EXPECT_EQ(0, rints[10]);
        EXPECT_EQ(1, rints[11]);
        EXPECT_EQ(0, rints[12]);
        EXPECT_EQ(1, rints[13]);

        EXPECT_SOFTEQ(rdbls[3], 0.63, 1.e-12);
        EXPECT_EQ(State::PLUS_X, rints[14]);
        EXPECT_EQ(1, rints[15]);

        EXPECT_SOFTEQ(rdbls[4], 0.63, 1.e-12);
        EXPECT_EQ(State::MINUS_X, rints[16]);
        EXPECT_EQ(1, rints[17]);
    }

    {
        single_shell_kernel2<<<1,1>>>(p2);
    }
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/RTK_Cell_Tester.cu
//---------------------------------------------------------------------------//
