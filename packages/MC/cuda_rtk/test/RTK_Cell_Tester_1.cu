//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_Cell_Tester_1.cu
 * \author Tom Evans
 * \date   Tue Nov 29 17:08:48 2016
 * \brief  RTK_Cell_Tester member and kernel definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/fill.h>

#include "CudaUtils/cuda_utils/Device_View.hh"

#include "RTK_Cell_Tester.hh"
#include "../RTK_Cell.cuh"
#include "../RTK_State_Vector.cuh"
#include "../RTK_State.cuh"

//---------------------------------------------------------------------------//
// TYPES
//---------------------------------------------------------------------------//

using Device_Cell   = cuda_profugus::RTK_Cell;
using Device_View   = cuda::Device_View<Device_Cell>;
using Manager       = cuda_profugus::RTK_Cell_DMM;
using Managers      = Device_View::Managers;
using DVF           = Device_View::DVF;
using Vector        = Device_Cell::Space_Vector;
using State         = cuda_profugus::RTK_State;
using State_Vec     = cuda_profugus::RTK_State_Vector;
using State_Vec_DMM = cuda_profugus::RTK_State_Vector_DMM;

using def::X; using def::Y; using def::Z;

//---------------------------------------------------------------------------//
// SINGLE-SHELL
//---------------------------------------------------------------------------//

__global__
void single_shell_kernel1(
    DVF          pins,
    State_Vec    states,
    int         *ints,
    double      *dbls,
    Vector      *svs)
{
    auto &pin = pins[0];

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

    int tid = cuda::utility::thread_id();

    // Make the state and initialize
    {
        pin.initialize({0.0, 0.53, 0.0}, states, tid);
        ints[6] = states.region(tid);
        ints[7] = states.segment(tid);
        pin.initialize({0.0, 0.54, 0.0}, states, tid);
        ints[8]  = states.region(tid);
        ints[9]  = states.segment(tid);
        ints[10] = pin.cell(states.region(tid), states.segment(tid));
        pin.initialize({0.0, 0.55, 0.0}, states, tid);
        ints[11] = states.region(tid);
        ints[12] = states.segment(tid);
        ints[13] = pin.cell(states.region(tid), states.segment(tid));
    }

    // Distance-to-boundary box-boundary tests
    Vector r, omega;
    {
        r     = {0.0, 0.59, 0.0};
        omega = {1.0, 0.0, 0.0};
        pin.initialize(r, states, tid);
        pin.distance_to_boundary(r, omega, states, tid);

        dbls[3]  = states.dist_to_next_region(tid);
        ints[14] = states.exiting_face(tid);
        ints[15] = states.region(tid);
    }
    {
        r     = {0.0, 0.59, 0.0};
        omega = {-1.0, 0.0, 0.0};
        pin.initialize(r, states, tid);
        pin.distance_to_boundary(r, omega, states, tid);

        dbls[4]  = states.dist_to_next_region(tid);
        ints[16] = states.exiting_face(tid);
        ints[17] = states.region(tid);
    }

    // Should not have vessel
    ints[18] = pin.has_vessel();
}

//---------------------------------------------------------------------------//

__global__
void single_shell_kernel2(
    Device_Cell  pin,
    State_Vec    states,
    int         *ints,
    double      *dbls)
{
    DEVICE_REQUIRE(ints[0] == -1);
    DEVICE_REQUIRE(dbls[0] == -1.0);

    Vector r, omega;
    int    n = 0, m = 0;

    int tid = cuda::utility::thread_id();

    // Pin intersection tests
    {
        r     = {  0.43,   0.51,   1.20};
        omega = { -0.07450781,  -0.17272265,   0.98214840};
        pin.initialize(r, states, tid);
        ints[n++] = states.region(tid);
        pin.distance_to_boundary(r, omega, states, tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.next_region(tid);
        dbls[m++] = states.dist_to_next_region(tid);
    }
    {
        r     = {  0.43,   0.51,   1.20};
        omega = {  0.01923789,  -0.98113214,  -0.19237885};
        pin.initialize(r, states, tid);
        ints[n++] = states.region(tid);
        pin.distance_to_boundary(r, omega, states, tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.next_region(tid);
        dbls[m++] = states.dist_to_next_region(tid);
    }
    {
        r     = { -0.49,   0.10,   1.20};
        omega = {  0.04377546,  -0.01122448,   0.99897834};
        pin.initialize(r, states, tid);
        ints[n++] = states.region(tid);
        pin.distance_to_boundary(r, omega, states, tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.next_region(tid);
        dbls[m++] = states.dist_to_next_region(tid);
    }
    {
        r     = { -0.31,  -0.46,   1.20};
        omega = { -0.01642565,   0.05812155,   0.99817438};
        pin.initialize(r, states, tid);
        ints[n++] = states.region(tid);
        pin.distance_to_boundary(r, omega, states, tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.next_region(tid);
        dbls[m++] = states.dist_to_next_region(tid);
    }
    {
        r     = {  0.21,  -0.56,   1.20};
        omega = { -0.03911234,   0.07065454,   0.99673375};
        pin.initialize(r, states, tid);
        ints[n++] = states.region(tid);
        pin.distance_to_boundary(r, omega, states, tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.next_region(tid);
        dbls[m++] = states.dist_to_next_region(tid);
    }
    {
        r     = {  0.21,  -0.56,  14.20};
        omega = {  0.021262916894,   0.051770580263,   0.998432619351};
        pin.initialize(r, states, tid);
        ints[n++] = states.region(tid);
        pin.distance_to_boundary(r, omega, states, tid);
        ints[n++] = states.exiting_face(tid);
        dbls[m++] = states.dist_to_next_region(tid);
    }

    {
        r     = {  0.42,   0.10,  14.20};
        omega = { -0.262232636986,  -0.030141682412,  -0.964533837188};
        pin.initialize(r, states, tid);
        ints[n++] = states.region(tid);
        pin.distance_to_boundary(r, omega, states, tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.next_region(tid);
        dbls[m++] = states.dist_to_next_region(tid);
    }
    {
        r     = {  0.10,   0.30,  12.10};
        omega = {  0.408248290464,   0.163299316186,  -0.898146239020};
        pin.initialize(r, states, tid);
        ints[n++] = states.region(tid);
        pin.distance_to_boundary(r, omega, states, tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.next_region(tid);
        dbls[m++] = states.dist_to_next_region(tid);
    }
    {
        r     = {  0.10,   0.30,  12.10};
        omega = {  0.011875513070,   0.004750205228,  -0.999918200524};
        pin.initialize(r, states, tid);
        ints[n++] = states.region(tid);
        pin.distance_to_boundary(r, omega, states, tid);
        ints[n++] = states.exiting_face(tid);
        dbls[m++] = states.dist_to_next_region(tid);
    }
    {
        r     = {  0.10,   0.30,  12.10};
        omega = {  0.072244077132,   0.028897630853,   0.996968264415};
        pin.initialize(r, states, tid);
        ints[n++] = states.region(tid);
        pin.distance_to_boundary(r, omega, states, tid);
        ints[n++] = states.exiting_face(tid);
        dbls[m++] = states.dist_to_next_region(tid);
    }
}

//---------------------------------------------------------------------------//

void Single_Shell::construct()
{
    Managers managers = {std::make_shared<Manager>(*pins[0]),
                         std::make_shared<Manager>(*pins[1])};
    Device_View dv(managers);

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);
    thrust::device_vector<Vector> svs(20);

    State_Vec_DMM states;
    states.initialize(1);

    single_shell_kernel1<<<1,1>>>(
        dv.get_view(), states.device_instance(), ints.data().get(),
        dbls.data().get(), svs.data().get());

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

    EXPECT_FALSE(rints[18]);
}

//---------------------------------------------------------------------------//

void Single_Shell::track()
{
    Manager dmm(*pins[1]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);

    State_Vec_DMM states;
    states.initialize(1);

    single_shell_kernel2<<<1,1>>>(
        pin, states.device_instance(), ints.data().get(), dbls.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    int    n   = 0, m = 0;
    double eps = 1.0e-6;

    EXPECT_EQ(1, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 1.2334036420, eps);
    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 0.41448110826, eps);
    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 1.1101358881, eps);
    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 3.4103552300, eps);
    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 1.860292469, eps);
    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 0.08 / 0.998432619351, eps);
    EXPECT_EQ(State::PLUS_Z, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 3.3176648414, eps);
    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(1, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 3.9914694397e-01, eps);
    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(1, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], -12.10 / -0.999918200524, eps);
    EXPECT_EQ(State::MINUS_Z, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_SOFTEQ(rdbls[m++], 2.18 / 0.996968264415, eps);
    EXPECT_EQ(State::PLUS_Z, rints[n++]);
}

//---------------------------------------------------------------------------//
// MULTI_SHELL
//---------------------------------------------------------------------------//

__global__
void multi_shell_kernel1(
    Device_Cell  pin,
    int         *ints,
    double      *dbls)
{
    int n = 0, m = 0;

    auto r = pin.radii();
    dbls[m++] = r[0];
    dbls[m++] = r[1];

    dbls[m++] = pin.pitch(0);
    dbls[m++] = pin.pitch(1);
    dbls[m++] = pin.height();

    ints[n++] = pin.num_regions();
    ints[n++] = pin.num_shells();
    ints[n++] = pin.num_segments();
    ints[n++] = pin.matid(0);
    ints[n++] = pin.matid(1);
    ints[n++] = pin.matid(2);

    ints[n++] = pin.region(0.1, 0.48);
    ints[n++] = pin.region(0.1, 0.479);
    ints[n++] = pin.region(0.5, 0.3);
    ints[n++] = pin.region(0.4, 0.35);
}

//---------------------------------------------------------------------------//

__global__
void multi_shell_kernel2(
    Device_Cell  pin,
    State_Vec    states,
    int         *ints,
    double      *dbls)
{
    int n = 0, m = 0;

    Vector r, omega;

    int tid = cuda::utility::thread_id();

    // Tracking tests
    r     = {  0.50,   0.30,  12.10};
    omega = { -0.740797197487,  -0.642024237822,   0.197545919330};
    pin.initialize(r, states, tid);
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.next_region(tid);

    pin.cross_surface(states, tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.region(tid);

    for (int i = 0; i < 3; ++i)
        r[i] += states.dist_to_next_region(tid) * omega[i];

    pin.distance_to_boundary(r, omega, states, tid);

    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.next_region(tid);

    pin.cross_surface(states, tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.region(tid);

    for (int i = 0; i < 3; ++i)
        r[i] += states.dist_to_next_region(tid) * omega[i];

    pin.distance_to_boundary(r, omega, states, tid);

    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.next_region(tid);

    pin.cross_surface(states, tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.region(tid);

    for (int i = 0; i < 3; ++i)
        r[i] += states.dist_to_next_region(tid) * omega[i];

    pin.distance_to_boundary(r, omega, states, tid);

    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.next_region(tid);

    pin.cross_surface(states, tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.region(tid);

    for (int i = 0; i < 3; ++i)
        r[i] += states.dist_to_next_region(tid) * omega[i];

    pin.distance_to_boundary(r, omega, states, tid);

    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.region(tid);
}

//---------------------------------------------------------------------------//

void Multi_Shell::construct()
{
    Manager dmm(*pins[0]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);

    multi_shell_kernel1<<<1,1>>>(
        pin, ints.data().get(), dbls.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    int n = 0, m = 0;

    EXPECT_SOFT_EQ(0.49, rdbls[m++]);
    EXPECT_SOFT_EQ(0.54, rdbls[m++]);

    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(1, rints[n++]);

    EXPECT_EQ(1.26,  rdbls[m++]);
    EXPECT_EQ(1.26,  rdbls[m++]);
    EXPECT_EQ(14.28, rdbls[m++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(3, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
}

//---------------------------------------------------------------------------//

void Multi_Shell::track()
{
    Manager dmm(*pins[0]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);

    State_Vec_DMM states;
    states.initialize(1);

    multi_shell_kernel2<<<1,1>>>(
        pin, states.device_instance(), ints.data().get(), dbls.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    int    n   = 0, m = 0;
    double eps = 1.0e-6;

    EXPECT_SOFTEQ(rdbls[m++], 0.0446878772402, eps);
    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(1, rints[n++]);

    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(1, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 0.0520128055639, eps);
    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 0.978336739656, eps);
    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(1, rints[n++]);

    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(1, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 0.0520128055639, eps);
    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);

    EXPECT_EQ(State::INTERNAL, rints[n++]);
    EXPECT_EQ(2, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 0.32149321506, eps);
    EXPECT_EQ(State::MINUS_Y, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
}

//---------------------------------------------------------------------------//

__global__
void multi_shell_kernel3(
    Device_Cell  pin,
    int         *ints,
    double      *dbls)
{
    int n = 0;

    ints[n++] = pin.num_cells();
    ints[n++] = pin.num_regions();
    ints[n++] = pin.num_shells();
    ints[n++] = pin.num_segments();

    ints[n++] = pin.region(0.28, 0.0);
    ints[n++] = pin.segment(0.28, 0.0);
    ints[n++] = pin.segment(0.28, -0.01);

    ints[n++] = pin.cell(0, 0);
    ints[n++] = pin.cell(1, 0);
    ints[n++] = pin.cell(2, 0);
    ints[n++] = pin.cell(3, 0);

    ints[n++] = pin.cell(0, 1);
    ints[n++] = pin.cell(1, 1);
    ints[n++] = pin.cell(2, 1);
    ints[n++] = pin.cell(3, 1);

    ints[n++] = pin.cell(0, 2);
    ints[n++] = pin.cell(1, 2);
    ints[n++] = pin.cell(2, 2);
    ints[n++] = pin.cell(3, 2);

    ints[n++] = pin.cell(0, 3);
    ints[n++] = pin.cell(1, 3);
    ints[n++] = pin.cell(2, 3);
    ints[n++] = pin.cell(3, 3);
}

//---------------------------------------------------------------------------//

__global__
void multi_shell_kernel4(
    Device_Cell  pin,
    State_Vec    states,
    int         *ints,
    double      *dbls)
{
    int n = 0, m = 0;

    Vector r, omega;

    int tid = cuda::utility::thread_id();

    // test tracking
    r     = {  0.43,   0.51,   1.20};
    omega = { -0.267261241912,  -0.534522483825,   0.801783725737};

    // initialize
    pin.initialize(r, states, tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.face(tid);

    // get distance to boundary
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
    ints[n++] = states.next_segment(tid);
    ints[n++] = states.exiting_face(tid);

    // move to surface
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    // cross the surface
    pin.cross_surface(states, tid);

    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.exiting_face(tid);

    // get distance to boundary
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
    ints[n++] = states.next_segment(tid);
    ints[n++] = states.exiting_face(tid);

    // move to surface
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    // cross the surface
    pin.cross_surface(states, tid);

    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.exiting_face(tid);

    // get distance to boundary
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
    ints[n++] = states.next_segment(tid);
    ints[n++] = states.exiting_face(tid);

    // move to surface
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    // cross the surface
    pin.cross_surface(states, tid);

    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.exiting_face(tid);

    // get distance to boundary
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
    ints[n++] = states.next_segment(tid);
    ints[n++] = states.exiting_face(tid);

    // move to surface
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    // cross the surface
    pin.cross_surface(states, tid);

    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.exiting_face(tid);

    // get distance to boundary
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
    ints[n++] = states.next_segment(tid);
    ints[n++] = states.exiting_face(tid);

    // move to surface
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    // cross the surface
    pin.cross_surface(states, tid);

    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.exiting_face(tid);

    // get distance to boundary
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
    ints[n++] = states.next_segment(tid);
    ints[n++] = states.exiting_face(tid);

    // move to surface
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    // cross the surface
    pin.cross_surface(states, tid);

    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.exiting_face(tid);

    // get distance to boundary
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
    ints[n++] = states.next_segment(tid);
    ints[n++] = states.exiting_face(tid);

    // move to surface
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    // cross the surface
    pin.cross_surface(states, tid);

    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.exiting_face(tid);

    // get distance to boundary
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
    ints[n++] = states.next_segment(tid);
    ints[n++] = states.exiting_face(tid);

    // move to surface
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    // cross the surface
    pin.cross_surface(states, tid);

    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.exiting_face(tid);

    // get distance to boundary
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.exiting_face(tid);
}

//---------------------------------------------------------------------------//

void Multi_Shell::multiseg_construct()
{
    Manager dmm(*pins[1]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);

    multi_shell_kernel3<<<1,1>>>(
        pin, ints.data().get(), dbls.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    int n = 0;

    EXPECT_EQ(16, rints[n++]);
    EXPECT_EQ(4,  rints[n++]);
    EXPECT_EQ(3,  rints[n++]);
    EXPECT_EQ(4,  rints[n++]);

    EXPECT_EQ(1,  rints[n++]);
    EXPECT_EQ(0,  rints[n++]);
    EXPECT_EQ(2,  rints[n++]);

    EXPECT_EQ(0,  rints[n++]);
    EXPECT_EQ(1,  rints[n++]);
    EXPECT_EQ(2,  rints[n++]);
    EXPECT_EQ(3,  rints[n++]);

    EXPECT_EQ(4,  rints[n++]);
    EXPECT_EQ(5,  rints[n++]);
    EXPECT_EQ(6,  rints[n++]);
    EXPECT_EQ(7,  rints[n++]);

    EXPECT_EQ(8,  rints[n++]);
    EXPECT_EQ(9,  rints[n++]);
    EXPECT_EQ(10, rints[n++]);
    EXPECT_EQ(11, rints[n++]);

    EXPECT_EQ(12, rints[n++]);
    EXPECT_EQ(13, rints[n++]);
    EXPECT_EQ(14, rints[n++]);
    EXPECT_EQ(15, rints[n++]);
}

//---------------------------------------------------------------------------//

void Multi_Shell::multiseg_track()
{
    Manager dmm(*pins[1]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(100, -1);
    thrust::device_vector<double> dbls(100, -1);

    State_Vec_DMM states;
    states.initialize(1);

    multi_shell_kernel4<<<1,1>>>(
        pin, states.device_instance(), ints.data().get(), dbls.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    int    n   = 0, m = 0;
    double eps = 1.0e-6;

    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(State::NONE, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 2.2028008712e-01, eps);
    EXPECT_EQ(State::NONE, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 9.4898743120e-02, eps);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 4.0177140025e-01, eps);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 2.3717240314e-01, eps);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(4, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_EQ(4, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 4.9908842021e-01, eps);
    EXPECT_EQ(4, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 1.5570162247e-01, eps);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 2.4606977777e-01, eps);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 9.4898743120e-02, eps);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(State::INTERNAL, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 1.8286351326e-01, eps);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(State::NONE, rints[n++]);
    EXPECT_EQ(State::MINUS_Y, rints[n++]);
}

//---------------------------------------------------------------------------//
// EMPTY CELL
//---------------------------------------------------------------------------//

__global__
void empty_kernel1(
    Device_Cell  box,
    State_Vec    states,
    int         *ints,
    double      *dbls)
{
    int n = 0, m = 0;

    ints[n++] = box.num_regions();
    ints[n++] = box.num_shells();

    ints[n++] = box.region(0.0, 0.5401);
    ints[n++] = box.region(0.0, 0.5399);

    ints[n++] = box.matid(0);

    int tid = cuda::utility::thread_id();

    Vector r, omega;
    {
        r     = {0.0, 0.59, 0.0};
        omega = {1.0, 0.0, 0.0};
        box.initialize(r, states, tid);
        box.distance_to_boundary(r, omega, states, tid);
        dbls[m++] = states.dist_to_next_region(tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.region(tid);
    }
    {
        r     = {0.0, 0.59, 0.0};
        omega = {-1.0, 0.0, 0.0};
        box.initialize(r, states, tid);
        box.distance_to_boundary(r, omega, states, tid);
        dbls[m++] = states.dist_to_next_region(tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.region(tid);
    }
    {
        r     = {  0.23,  -0.63,  12.10};
        omega = { -0.591113929288,   0.783346101414,   0.192232172126};
        box.initialize(r, states, tid);
        box.distance_to_boundary(r, omega, states, tid);
        dbls[m++] = states.dist_to_next_region(tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.region(tid);
    }
    {
        r     = {  0.45,  -0.62,  12.10};
        omega = { -0.628969195431,   0.754763034517,   0.186361243091};
        box.initialize(r, states, tid);
        box.distance_to_boundary(r, omega, states, tid);
        dbls[m++] = states.dist_to_next_region(tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.region(tid);
    }
    {
        r     = { -0.12,  -0.38,  12.10};
        omega = {  0.810238620974,  -0.492497985298,   0.317740635676};
        box.initialize(r, states, tid);
        box.distance_to_boundary(r, omega, states, tid);
        dbls[m++] = states.dist_to_next_region(tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.region(tid);
    }
    {
        r     = { -0.12,  -0.38,  12.10};
        omega = {  0.097437248619,  -0.059226562886,   0.993477829058};
        box.initialize(r, states, tid);
        box.distance_to_boundary(r, omega, states, tid);
        dbls[m++] = states.dist_to_next_region(tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.region(tid);
    }
    {
        r     = { -0.12,  -0.38,   2.10};
        omega = {  0.080591552144,  -0.048987021891,  -0.995542702956};
        box.initialize(r, states, tid);
        box.distance_to_boundary(r, omega, states, tid);
        dbls[m++] = states.dist_to_next_region(tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.region(tid);
    }
}

//---------------------------------------------------------------------------//

__global__
void empty_kernel2(
    Device_Cell  box,
    State_Vec    states,
    int         *ints,
    double      *dbls)
{
    int n = 0, m = 0;

    ints[n++] = box.num_regions();
    ints[n++] = box.num_shells();

    ints[n++] = box.region(0.0, 0.5401);
    ints[n++] = box.region(0.0, 0.5399);

    ints[n++] = box.matid(0);

    dbls[m++] = box.pitch(0);
    dbls[m++] = box.pitch(1);

    int tid = cuda::utility::thread_id();

    Vector r, omega;
    {
        r     = {0.0, 0.59, 0.0};
        omega = {1.0, 0.0, 0.0};
        box.initialize(r, states, tid);
        box.distance_to_boundary(r, omega, states, tid);
        dbls[m++] = states.dist_to_next_region(tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.region(tid);
    }
    {
        r     = {0.0, 0.59, 0.0};
        omega = {-1.0, 0.0, 0.0};
        box.initialize(r, states, tid);
        box.distance_to_boundary(r, omega, states, tid);
        dbls[m++] = states.dist_to_next_region(tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.region(tid);
    }
    {
        r     = {  0.123,  -0.63,  12.10};
        omega = { -0.556102184819,   0.807165237093,   0.198077358796};
        box.initialize(r, states, tid);
        box.distance_to_boundary(r, omega, states, tid);
        dbls[m++] = states.dist_to_next_region(tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.region(tid);
    }
    {
        r     = {  0.045,  -0.62,  12.10};
        omega = { -0.080640139920,   0.967681679037,   0.238933747910};
        box.initialize(r, states, tid);
        box.distance_to_boundary(r, omega, states, tid);
        dbls[m++] = states.dist_to_next_region(tid);
        ints[n++] = states.exiting_face(tid);
        ints[n++] = states.region(tid);
    }
}

//---------------------------------------------------------------------------//

void Empty::square()
{
    Manager dmm(*pins[0]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);

    State_Vec_DMM states;
    states.initialize(1);

    empty_kernel1<<<1,1>>>(
        pin, states.device_instance(), ints.data().get(), dbls.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    int    n   = 0, m = 0;
    double eps = 1.0e-6;

    EXPECT_EQ(1, ints[n++]);
    EXPECT_EQ(0, ints[n++]);

    EXPECT_EQ(0, ints[n++]);
    EXPECT_EQ(0, ints[n++]);

    EXPECT_EQ(11, ints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 0.63, 1.e-12);
    EXPECT_EQ(State::PLUS_X, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 0.63, 1.e-12);
    EXPECT_EQ(State::MINUS_X, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 1.4548802818, eps);
    EXPECT_EQ(State::MINUS_X, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 1.6561489406, eps);
    EXPECT_EQ(State::PLUS_Y, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 0.50761628974, eps);
    EXPECT_EQ(State::MINUS_Y, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 2.1943116758, eps);
    EXPECT_EQ(State::PLUS_Z, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 2.1094022323, eps);
    EXPECT_EQ(State::MINUS_Z, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
}

//---------------------------------------------------------------------------//

void Empty::rectangle()
{
    Manager dmm(*pins[1]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);

    State_Vec_DMM states;
    states.initialize(1);

    empty_kernel2<<<1,1>>>(
        pin, states.device_instance(), ints.data().get(), dbls.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    int    n   = 0, m = 0;
    double eps = 1.0e-6;

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_EQ(11, rints[n++]);

    EXPECT_EQ(0.25, rdbls[m++]);
    EXPECT_EQ(1.26, rdbls[m++]);

    EXPECT_SOFTEQ(rdbls[m++], 0.125, 1.e-12);
    EXPECT_EQ(State::PLUS_X, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 0.125, 1.e-12);
    EXPECT_EQ(State::MINUS_X, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 4.4596120420e-01, eps);
    EXPECT_EQ(State::MINUS_X, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_SOFTEQ(rdbls[m++], 1.2917470973e+00, eps);
    EXPECT_EQ(State::PLUS_Y, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
}

//---------------------------------------------------------------------------//
// end of RTK_Cell_Tester_1.cu
//---------------------------------------------------------------------------//
