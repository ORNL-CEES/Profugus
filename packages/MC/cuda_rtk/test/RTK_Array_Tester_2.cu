//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_Array_Tester_2.cu
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

    ints[m++] = state.level_coord[0][0];
    ints[m++] = state.level_coord[0][1];
    ints[m++] = state.level_coord[0][2];
    ints[m++] = state.region;
    ints[m++] = array.matid(state);

    r = {1.259, 1.27, 1.1};
    array.initialize(r, state);

    ints[m++] = state.level_coord[0][0];
    ints[m++] = state.level_coord[0][1];
    ints[m++] = state.level_coord[0][2];
    ints[m++] = state.region;
    ints[m++] = array.matid(state);

    r = {3.560000,   2.239887,   1.300000};
    array.initialize(r, state);

    ints[m++] = state.level_coord[0][0];
    ints[m++] = state.level_coord[0][1];
    ints[m++] = state.level_coord[0][2];
    ints[m++] = state.region;
    ints[m++] = array.matid(state);

    r = {1.570000,   0.931993,   2.700000};
    array.initialize(r, state);

    ints[m++] = state.level_coord[0][0];
    ints[m++] = state.level_coord[0][1];
    ints[m++] = state.level_coord[0][2];
    ints[m++] = state.region;
    ints[m++] = array.matid(state);

    r = {1.300000,   2.044919,   3.800000};
    array.initialize(r, state);

    ints[m++] = state.level_coord[0][0];
    ints[m++] = state.level_coord[0][1];
    ints[m++] = state.level_coord[0][2];
    ints[m++] = state.region;
    ints[m++] = array.matid(state);
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

    EXPECT_EQ(1,  rints[m++]);
    EXPECT_EQ(1,  rints[m++]);
    EXPECT_EQ(0,  rints[m++]);
    EXPECT_EQ(1,  rints[m++]);
    EXPECT_EQ(10, rints[m++]);

    EXPECT_EQ(0,  rints[m++]);
    EXPECT_EQ(1,  rints[m++]);
    EXPECT_EQ(0,  rints[m++]);
    EXPECT_EQ(1,  rints[m++]);
    EXPECT_EQ(10, rints[m++]);

    EXPECT_EQ(2,  rints[m++]);
    EXPECT_EQ(1,  rints[m++]);
    EXPECT_EQ(0,  rints[m++]);
    EXPECT_EQ(0,  rints[m++]);
    EXPECT_EQ(5,  rints[m++]);

    EXPECT_EQ(1,  rints[m++]);
    EXPECT_EQ(0,  rints[m++]);
    EXPECT_EQ(0,  rints[m++]);
    EXPECT_EQ(0,  rints[m++]);
    EXPECT_EQ(3,  rints[m++]);

    EXPECT_EQ(1,  rints[m++]);
    EXPECT_EQ(1,  rints[m++]);
    EXPECT_EQ(0,  rints[m++]);
    EXPECT_EQ(0,  rints[m++]);
    EXPECT_EQ(0,  rints[m++]);
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

    ints[m++] = state.level_coord[1][0];
    ints[m++] = state.level_coord[1][1];
    ints[m++] = state.level_coord[1][2];
    ints[m++] = state.level_coord[0][0];
    ints[m++] = state.level_coord[0][1];
    ints[m++] = state.level_coord[0][2];
    ints[m++] = state.region;

    r = {8.34, 2.3, -3.1};
    array.initialize(r, state);

    ints[m++] = state.level_coord[1][0];
    ints[m++] = state.level_coord[1][1];
    ints[m++] = state.level_coord[1][2];
    ints[m++] = state.level_coord[0][0];
    ints[m++] = state.level_coord[0][1];
    ints[m++] = state.level_coord[0][2];
    ints[m++] = state.region;
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

    EXPECT_EQ(0, rints[m++]);
    EXPECT_EQ(0, rints[m++]);
    EXPECT_EQ(2, rints[m++]);
    EXPECT_EQ(0, rints[m++]);
    EXPECT_EQ(0, rints[m++]);
    EXPECT_EQ(0, rints[m++]);
    EXPECT_EQ(1, rints[m++]);

    EXPECT_EQ(2, rints[m++]);
    EXPECT_EQ(0, rints[m++]);
    EXPECT_EQ(0, rints[m++]);
    EXPECT_EQ(0, rints[m++]);
    EXPECT_EQ(1, rints[m++]);
    EXPECT_EQ(0, rints[m++]);
    EXPECT_EQ(0, rints[m++]);
}

//---------------------------------------------------------------------------//
// CORE
//---------------------------------------------------------------------------//

__global__
void regcoreA_kernel(
    Core_Array  core,
    int        *ints,
    double     *dbls)
{
    using def::X; using def::Y; using def::Z;

    State state;
    Vector r, omega;
    int m = 0, n = 0;
    double d = 0.0;

    r     = {  0.044859500000,   5.638180000000,   7.185140000000};
    omega = {  0.994391000000,   0.099447500000,  -0.036019600000};
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

    core.cross_surface(r, state);

    ints[n++] = state.level_coord[1][X];
    ints[n++] = state.level_coord[1][Y];
    ints[n++] = state.level_coord[1][Z];
    ints[n++] = state.level_coord[0][X];
    ints[n++] = state.level_coord[0][Y];
    ints[n++] = state.level_coord[0][Z];
    ints[n++] = state.escaping_face;
}

//---------------------------------------------------------------------------//

__global__
void regcoreB_kernel(
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

    // escape
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
    ints[n++] = state.escaping_face;
}

//---------------------------------------------------------------------------//

void RegCore::run_test()
{
    // Make DMM
    Core_Manager dmm(*core);

    // Get the host object
    auto array = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);

    regcoreA_kernel<<<1,1>>>(array, ints.data().get(), dbls.data().get());

    int m = 0, n = 0;
    double eps = 1.0e-6;

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(5, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 2.489101872, eps);
    EXPECT_EQ(State::PLUS_X, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(5, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 2.534214409, eps);
    EXPECT_EQ(State::PLUS_X, rints[n++]);

    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(5, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 2.534214409, eps);
    EXPECT_EQ(State::PLUS_X, rints[n++]);

    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(State::PLUS_X, rints[n++]);

    regcoreB_kernel<<<1,1>>>(array, ints.data().get(), dbls.data().get());

    rints = thrust::host_vector<int>(ints.begin(), ints.end());
    rdbls = thrust::host_vector<double>(dbls.begin(), dbls.end());

    m = 0;
    n = 0;

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

    EXPECT_EQ(2,  rints[n++]);
    EXPECT_EQ(2,  rints[n++]);
    EXPECT_EQ(-1, rints[n++]);
    EXPECT_EQ(0,  rints[n++]);
    EXPECT_EQ(0,  rints[n++]);
    EXPECT_EQ(-1, rints[n++]);
    EXPECT_EQ(State::MINUS_Z, rints[n++]);
}

//---------------------------------------------------------------------------//
// end of RTK_Array_Tester_2.cu
//---------------------------------------------------------------------------//
