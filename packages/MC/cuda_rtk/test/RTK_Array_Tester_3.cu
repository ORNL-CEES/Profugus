//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_Array_Tester_3.cu
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
// CORE
//---------------------------------------------------------------------------//

__global__
void reflect_kernel(
    Core_Array  core,
    int        *ints,
    double     *dbls)
{
    using def::X; using def::Y; using def::Z;

    State state;
    Vector r, omega;
    int m = 0, n = 0;
    double d = 0.0;

    r     = {  4.202350000000,   2.820900000000,  18.507800000000};
    omega = {  0.098705500000,   0.137387000000,  -0.985587000000};

    core.initialize(r, state);
    core.distance_to_boundary(r, omega, state);

    d = state.dist_to_next_region;

    ints[n++] = state.level_coord[1][X];
    ints[n++] = state.level_coord[1][Y];
    ints[n++] = state.level_coord[1][Z];
    ints[n++] = state.level_coord[0][X];
    ints[n++] = state.level_coord[0][Y];
    ints[n++] = state.level_coord[0][Z];

    ints[n++] = state.region;
    ints[n++] = core.matid(state);
    dbls[m++] = d;
    ints[n++] = state.exiting_face;

    // process surface
    r[0] = r[0] + d * omega[0];
    r[1] = r[1] + d * omega[1];
    r[2] = r[2] + d * omega[2];

    core.cross_surface(r, state);

    // next step
    core.distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;

    ints[n++] = state.level_coord[1][X];
    ints[n++] = state.level_coord[1][Y];
    ints[n++] = state.level_coord[1][Z];
    ints[n++] = state.level_coord[0][X];
    ints[n++] = state.level_coord[0][Y];
    ints[n++] = state.level_coord[0][Z];

    ints[n++] = state.region;
    ints[n++] = core.matid(state);
    dbls[m++] = d;
    ints[n++] = state.exiting_face;
    ints[n++] = state.next_region;
    ints[n++] = state.next_face;

    // process surface
    r[0] = r[0] + d * omega[0];
    r[1] = r[1] + d * omega[1];
    r[2] = r[2] + d * omega[2];

    core.cross_surface(r, state);
    ints[n++] = state.face;

    // next step
    core.distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;

    ints[n++] = state.level_coord[1][X];
    ints[n++] = state.level_coord[1][Y];
    ints[n++] = state.level_coord[1][Z];
    ints[n++] = state.level_coord[0][X];
    ints[n++] = state.level_coord[0][Y];
    ints[n++] = state.level_coord[0][Z];

    ints[n++] = state.region;
    ints[n++] = core.matid(state);
    dbls[m++] = d;
    ints[n++] = state.exiting_face;

    // process surface
    r[0] = r[0] + d * omega[0];
    r[1] = r[1] + d * omega[1];
    r[2] = r[2] + d * omega[2];

    core.cross_surface(r, state);

    // next step
    core.distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;

    ints[n++] = state.level_coord[1][X];
    ints[n++] = state.level_coord[1][Y];
    ints[n++] = state.level_coord[1][Z];
    ints[n++] = state.level_coord[0][X];
    ints[n++] = state.level_coord[0][Y];
    ints[n++] = state.level_coord[0][Z];

    ints[n++] = state.region;
    ints[n++] = core.matid(state);
    dbls[m++] = d;
    ints[n++] = state.exiting_face;

    // process surface
    r[0] = r[0] + d * omega[0];
    r[1] = r[1] + d * omega[1];
    r[2] = r[2] + d * omega[2];

    core.cross_surface(r, state);

    // next step
    core.distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;

    ints[n++] = state.level_coord[1][X];
    ints[n++] = state.level_coord[1][Y];
    ints[n++] = state.level_coord[1][Z];
    ints[n++] = state.level_coord[0][X];
    ints[n++] = state.level_coord[0][Y];
    ints[n++] = state.level_coord[0][Z];

    ints[n++] = state.region;
    ints[n++] = core.matid(state);
    dbls[m++] = d;
    ints[n++] = state.exiting_face;

    // process surface
    r[0] = r[0] + d * omega[0];
    r[1] = r[1] + d * omega[1];
    r[2] = r[2] + d * omega[2];

    core.cross_surface(r, state);

    // next step
    core.distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;

    ints[n++] = state.level_coord[1][X];
    ints[n++] = state.level_coord[1][Y];
    ints[n++] = state.level_coord[1][Z];
    ints[n++] = state.level_coord[0][X];
    ints[n++] = state.level_coord[0][Y];
    ints[n++] = state.level_coord[0][Z];

    ints[n++] = state.region;
    ints[n++] = core.matid(state);
    dbls[m++] = d;
    ints[n++] = state.exiting_face;

    // reflecting
    r[0] = r[0] + d * omega[0];
    r[1] = r[1] + d * omega[1];
    r[2] = r[2] + d * omega[2];

    core.cross_surface(r, state);

    ints[n++] = state.level_coord[1][X];
    ints[n++] = state.level_coord[1][Y];
    ints[n++] = state.level_coord[1][Z];
    ints[n++] = state.level_coord[0][X];
    ints[n++] = state.level_coord[0][Y];
    ints[n++] = state.level_coord[0][Z];
    ints[n++] = state.reflecting_face;
    ints[n++] = state.exiting_face;

    // flip direction and make sure that next-distance-to-boundary
    // clears the reflecting flag
    omega[2] = omega[2] - 2.0 * omega[2];
    core.distance_to_boundary(r, omega, state);
    d = state.dist_to_next_region;

    ints[n++] = state.level_coord[1][X];
    ints[n++] = state.level_coord[1][Y];
    ints[n++] = state.level_coord[1][Z];
    ints[n++] = state.level_coord[0][X];
    ints[n++] = state.level_coord[0][Y];
    ints[n++] = state.level_coord[0][Z];

    ints[n++] = state.region;
    ints[n++] = core.matid(state);
    dbls[m++] = d;
    ints[n++] = state.reflecting_face;
    ints[n++] = state.exiting_face;
}

//---------------------------------------------------------------------------//

void RegCore::reflect_test()
{
    // Make DMM
    Core_Manager dmm(*core);

    // Get the host object
    auto array = dmm.device_instance();

    thrust::device_vector<int>    ints(75, -1);
    thrust::device_vector<double> dbls(50, -1);

    reflect_kernel<<<1,1>>>(array, ints.data().get(), dbls.data().get());

    int m = 0, n = 0;
    double eps = 1.0e-6;

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(5, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 4.289626385, eps);
    EXPECT_EQ(State::MINUS_Z, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 1.195583643, eps);
    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(5, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 1.495799820, eps);
    EXPECT_EQ(State::PLUS_Y, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(5, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 1.505346029, eps);
    EXPECT_EQ(State::PLUS_X, rints[n++]);

    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(5, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 7.665827372, eps);
    EXPECT_EQ(State::PLUS_Y, rints[n++]);

    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(5, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 2.626270607, eps);
    EXPECT_EQ(State::MINUS_Z, rints[n++]);

    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(State::MINUS_Z, rints[n++]);
    EXPECT_EQ(State::MINUS_Z, rints[n++]);

    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(5, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 14.488827470, eps);
    EXPECT_EQ(State::NONE, rints[n++]);
    EXPECT_EQ(State::PLUS_Z, rints[n++]);
}

//---------------------------------------------------------------------------//
// end of RTK_Array_Tester_3.cu
//---------------------------------------------------------------------------//
