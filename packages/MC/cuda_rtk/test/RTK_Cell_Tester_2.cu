//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_Cell_Tester_2.cu
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
#include "../RTK_State.cuh"
#include "../RTK_State_Vector.cuh"

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
// GAP
//---------------------------------------------------------------------------//

__global__
void gap_kernel1(
    Device_Cell  pin,
    State_Vec    states,
    int         *ints,
    double      *dbls,
    Vector      *svs)
{
    int n = 0, m = 0;

    Vector r, omega;

    int tid = cuda::utility::thread_id();

    pin.get_extents(svs[0], svs[1]);

    dbls[m++] = pin.radii()[0];
    dbls[m++] = pin.radii()[1];
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
    ints[n++] = pin.region(-0.7, -0.7);
    ints[n++] = pin.region(-0.64, 0.1);
    ints[n++] = pin.region(0.1, -0.64);

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

__global__
void gap_kernel2(
    Device_Cell  pin,
    int         *ints,
    double      *dbls,
    Vector      *svs)
{
    int n = 0, m = 0;

    pin.get_extents(svs[0], svs[1]);

    dbls[m++] = pin.radii()[0];
    dbls[m++] = pin.radii()[1];
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
    ints[n++] = pin.region(-0.7, 0.7);
    ints[n++] = pin.region(-0.64, 0.1);
    ints[n++] = pin.region(0.1, 0.64);
}

//---------------------------------------------------------------------------//

__global__
void gap_kernel3(
    Device_Cell  pin,
    int         *ints,
    double      *dbls,
    Vector      *svs)
{
    int n = 0, m = 0;

    pin.get_extents(svs[0], svs[1]);

    dbls[m++] = pin.radii()[0];
    dbls[m++] = pin.radii()[1];
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
    ints[n++] = pin.region(0.7, 0.7);
    ints[n++] = pin.region(0.64, 0.1);
    ints[n++] = pin.region(0.1, 0.64);
}

//---------------------------------------------------------------------------//

__global__
void gap_kernel4(
    Device_Cell  pin,
    int         *ints,
    double      *dbls,
    Vector      *svs)
{
    int n = 0, m = 0;

    pin.get_extents(svs[0], svs[1]);

    dbls[m++] = pin.pitch(0);
    dbls[m++] = pin.pitch(1);
    dbls[m++] = pin.height();

    ints[n++] = pin.num_regions();
    ints[n++] = pin.num_shells();
    ints[n++] = pin.num_segments();

    ints[n++] = pin.matid(0);

    ints[n++] = pin.region(0.1, 0.48);
    ints[n++] = pin.region(0.1, 0.479);
    ints[n++] = pin.region(0.5, 0.3);
    ints[n++] = pin.region(0.4, 0.35);
    ints[n++] = pin.region(0.7, 0.7);
    ints[n++] = pin.region(0.64, 0.1);
    ints[n++] = pin.region(0.1, 0.64);

    ints[n++] = pin.has_vessel();
}

//---------------------------------------------------------------------------//

void Gap::lox_loy()
{
    Manager dmm(*pins[0]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);
    thrust::device_vector<Vector> svs(20);

    State_Vec_DMM states;
    states.initialize(1);

    gap_kernel1<<<1,1>>>(
        pin, states.device_instance(), ints.data().get(), dbls.data().get(),
        svs.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());
    thrust::host_vector<Vector> rsvs(svs.begin(), svs.end());

    int n = 0, m = 0;
    double eps = 1.0e-6;

    EXPECT_SOFT_EQ(-0.73, rsvs[0][0]);
    EXPECT_SOFT_EQ(-0.73, rsvs[0][1]);
    EXPECT_SOFT_EQ(0.0,   rsvs[0][2]);
    EXPECT_SOFT_EQ(0.63,  rsvs[1][0]);
    EXPECT_SOFT_EQ(0.63,  rsvs[1][1]);
    EXPECT_SOFT_EQ(14.28, rsvs[1][2]);

    EXPECT_SOFT_EQ(0.49, rdbls[m++]);
    EXPECT_SOFT_EQ(0.54, rdbls[m++]);
    EXPECT_SOFT_EQ(1.36, rdbls[m++]);
    EXPECT_SOFT_EQ(1.36, rdbls[m++]);
    EXPECT_EQ(14.28, rdbls[m++]);

    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(1, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(3, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);

    // Tracking results
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

    EXPECT_SOFTEQ(rdbls[m++], 0.4772505785762855, eps);
    EXPECT_EQ(State::MINUS_Y, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
}

//---------------------------------------------------------------------------//

void Gap::lox_hiy()
{
    Manager dmm(*pins[1]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);
    thrust::device_vector<Vector> svs(20);

    gap_kernel2<<<1,1>>>(
        pin, ints.data().get(), dbls.data().get(), svs.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());
    thrust::host_vector<Vector> rsvs(svs.begin(), svs.end());

    int n = 0, m = 0;

    EXPECT_SOFT_EQ(-0.73, rsvs[0][0]);
    EXPECT_SOFT_EQ(-0.63, rsvs[0][1]);
    EXPECT_SOFT_EQ(0.0,   rsvs[0][2]);
    EXPECT_SOFT_EQ(0.63,  rsvs[1][0]);
    EXPECT_SOFT_EQ(0.73,  rsvs[1][1]);
    EXPECT_SOFT_EQ(14.28, rsvs[1][2]);

    EXPECT_SOFT_EQ(0.49, rdbls[m++]);
    EXPECT_SOFT_EQ(0.54, rdbls[m++]);
    EXPECT_SOFT_EQ(1.36, rdbls[m++]);
    EXPECT_SOFT_EQ(1.36, rdbls[m++]);
    EXPECT_EQ(14.28, rdbls[m++]);

    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(1, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(3, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
}

//---------------------------------------------------------------------------//

void Gap::hix_hiy()
{
    Manager dmm(*pins[2]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);
    thrust::device_vector<Vector> svs(20);

    gap_kernel3<<<1,1>>>(
        pin, ints.data().get(), dbls.data().get(), svs.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());
    thrust::host_vector<Vector> rsvs(svs.begin(), svs.end());

    int n = 0, m = 0;

    EXPECT_SOFT_EQ(-0.63, rsvs[0][0]);
    EXPECT_SOFT_EQ(-0.63, rsvs[0][1]);
    EXPECT_SOFT_EQ(0.0,   rsvs[0][2]);
    EXPECT_SOFT_EQ(0.73,  rsvs[1][0]);
    EXPECT_SOFT_EQ(0.73,  rsvs[1][1]);
    EXPECT_SOFT_EQ(14.28, rsvs[1][2]);

    EXPECT_SOFT_EQ(0.49, rdbls[m++]);
    EXPECT_SOFT_EQ(0.54, rdbls[m++]);
    EXPECT_SOFT_EQ(1.36, rdbls[m++]);
    EXPECT_SOFT_EQ(1.36, rdbls[m++]);
    EXPECT_EQ(14.28, rdbls[m++]);

    EXPECT_EQ(3, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(1, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(3, rints[n++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
    EXPECT_EQ(2, rints[n++]);
}

//---------------------------------------------------------------------------//

void Gap::hom_hix_hiy()
{
    Manager dmm(*pins[3]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);
    thrust::device_vector<Vector> svs(20);

    gap_kernel4<<<1,1>>>(
        pin, ints.data().get(), dbls.data().get(), svs.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());
    thrust::host_vector<Vector> rsvs(svs.begin(), svs.end());

    int n = 0, m = 0;

    EXPECT_SOFT_EQ(-0.63, rsvs[0][0]);
    EXPECT_SOFT_EQ(-0.63, rsvs[0][1]);
    EXPECT_SOFT_EQ(0.0,   rsvs[0][2]);
    EXPECT_SOFT_EQ(0.73,  rsvs[1][0]);
    EXPECT_SOFT_EQ(0.73,  rsvs[1][1]);
    EXPECT_SOFT_EQ(14.28, rsvs[1][2]);

    EXPECT_SOFT_EQ(1.36, rdbls[m++]);
    EXPECT_SOFT_EQ(1.36, rdbls[m++]);
    EXPECT_EQ(14.28, rdbls[m++]);

    EXPECT_EQ(1, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(1, rints[n++]);

    EXPECT_EQ(3, rints[n++]);

    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);
    EXPECT_EQ(0, rints[n++]);

    EXPECT_FALSE(rints[n++]);
}

//---------------------------------------------------------------------------//
// VESSEL
//---------------------------------------------------------------------------//

__global__
void vessel_kernel1(
    Device_Cell  pin,
    State_Vec    states,
    int         *ints,
    double      *dbls)
{
    int n = 0, m = 0;

    int tid = cuda::utility::thread_id();

    Vector r, omega;

    ints[n++] = pin.has_vessel();

    r     = { -0.52,  -0.25,  11.80};
    omega = {  0.942738909161,   0.021633902593,  -0.332829270666};
    pin.initialize(r, states, tid);
    pin.distance_to_boundary(r, omega, states, tid);
    ints[n++] = states.region(tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);

    r     = { -2.52,  -4.20,  11.80};
    omega = {  0.130844809690,  -0.929686805689,  -0.344328446552};
    pin.initialize(r, states, tid);
    pin.distance_to_boundary(r, omega, states, tid);
    ints[n++] = states.region(tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);

    r     = { -0.52,  -0.25,  11.80};
    omega = { -0.638224253549,  -0.728778909543,  -0.248094947929};
    pin.initialize(r, states, tid);
    ints[n++] = states.exiting_face(tid);
    pin.distance_to_boundary(r, omega, states, tid);
    ints[n++] = states.region(tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);

    r     = { -2.52,  -4.20,  11.80};
    omega = {  0.077270894187,   0.976053400258,  -0.203344458387};
    pin.initialize(r, states, tid);
    pin.distance_to_boundary(r, omega, states, tid);
    ints[n++] = states.region(tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
}

//---------------------------------------------------------------------------//

__global__
void vessel_kernel2(
    Device_Cell  pin,
    State_Vec    states,
    int         *ints,
    double      *dbls)
{
    int n = 0, m = 0;

    int tid = cuda::utility::thread_id();

    Vector r, omega;

    ints[n++] = pin.has_vessel();

    r     = { -0.52,  -0.25,  11.80};
    omega = {  0.942738909161,   0.021633902593,  -0.332829270666};
    pin.initialize(r, states, tid);
    pin.distance_to_boundary(r, omega, states, tid);
    ints[n++] = states.region(tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);

    r     = { -2.52,  -4.20,  11.80};
    omega = {  0.130844809690,  -0.929686805689,  -0.344328446552};
    pin.initialize(r, states, tid);
    pin.distance_to_boundary(r, omega, states, tid);
    ints[n++] = states.region(tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);

    r     = { -0.52,  -0.25,  11.80};
    omega = { -0.638224253549,  -0.728778909543,  -0.248094947929};
    pin.initialize(r, states, tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.exiting_face(tid);
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);

    r     = { -2.52,  -4.20,  11.80};
    omega = {  0.077270894187,   0.976053400258,  -0.203344458387};
    pin.initialize(r, states, tid);
    pin.distance_to_boundary(r, omega, states, tid);
    ints[n++] = states.region(tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
}

//---------------------------------------------------------------------------//

__global__
void vessel_kernel3(
    Device_Cell  pin,
    State_Vec    states,
    int         *ints,
    double      *dbls)
{
    int n = 0, m = 0;

    int tid = cuda::utility::thread_id();

    Vector r, omega;

    ints[n++] = pin.has_vessel();

    r     = { -0.52,  -0.25,  11.80};
    omega = {  0.942738909161,   0.021633902593,  -0.332829270666};
    pin.initialize(r, states, tid);
    pin.distance_to_boundary(r, omega, states, tid);
    ints[n++] = states.region(tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);

    r     = { -0.52,  -0.25,  11.80};
    omega = { -0.638224253549,  -0.728778909543,  -0.248094947929};
    pin.initialize(r, states, tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.exiting_face(tid);
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);

    r     = { -2.52,  -4.20,  11.80};
    omega = {  0.130844809690,  -0.929686805689,  -0.344328446552};
    pin.initialize(r, states, tid);
    ints[n++] = states.region(tid);
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);

    r     = { -2.52,  -4.20,  11.80};
    omega = {  0.077270894187,   0.976053400258,  -0.203344458387};
    pin.initialize(r, states, tid);
    ints[n++] = states.region(tid);
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
}

//---------------------------------------------------------------------------//

__global__
void vessel_kernel4(
    Device_Cell  pin,
    State_Vec    states,
    int         *ints,
    double      *dbls)
{
    int n = 0, m = 0;

    int tid = cuda::utility::thread_id();

    Vector r, omega;

    ints[n++] = pin.has_vessel();

    r     = { -0.52,  -0.25,  11.80};
    omega = { -0.638224253549,  -0.728778909543,  -0.248094947929};
    pin.initialize(r, states, tid);
    ints[n++] = states.exiting_face(tid);
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
    ints[n++] = states.region(tid);
    ints[n++] = pin.matid(states.region(tid));

    // move the ray
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    // cross the surface
    pin.cross_surface(states, tid);

    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.exiting_face(tid);

    // next segment (to lower next shell)
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
    ints[n++] = states.region(tid);
    ints[n++] = pin.matid(states.region(tid));

    // move the ray
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    // cross the surface
    pin.cross_surface(states, tid);

    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.exiting_face(tid);

    // next segment (to low x-face)
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = pin.matid(states.region(tid));

    // move the ray
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    dbls[m++] = r[0];
}

//---------------------------------------------------------------------------//

__global__
void vessel_kernel5(
    Device_Cell  pin,
    State_Vec    states,
    int         *ints,
    double      *dbls)
{
    int n = 0, m = 0;

    int tid = cuda::utility::thread_id();

    Vector r, omega;

    ints[n++] = pin.has_vessel();

    r     = {-10.00,   0.25,  11.80};
    omega = {0.723129219960,   0.666634749651,  -0.180782304990};
    pin.initialize(r, states, tid);
    ints[n++] = states.exiting_face(tid);
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
    ints[n++] = states.region(tid);
    ints[n++] = pin.matid(states.region(tid));

    // move the ray
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    // cross the surface
    pin.cross_surface(states, tid);

    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.exiting_face(tid);

    // next segment (to lower next shell)
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.next_face(tid);
    ints[n++] = states.next_region(tid);
    ints[n++] = states.region(tid);
    ints[n++] = pin.matid(states.region(tid));

    // move the ray
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    // cross the surface
    pin.cross_surface(states, tid);

    ints[n++] = states.face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = states.segment(tid);
    ints[n++] = states.exiting_face(tid);

    // next segment (to low x-face)
    pin.distance_to_boundary(r, omega, states, tid);
    dbls[m++] = states.dist_to_next_region(tid);
    ints[n++] = states.exiting_face(tid);
    ints[n++] = states.region(tid);
    ints[n++] = pin.matid(states.region(tid));

    // move the ray
    r[0] += omega[0] * states.dist_to_next_region(tid);
    r[1] += omega[1] * states.dist_to_next_region(tid);
    r[2] += omega[2] * states.dist_to_next_region(tid);

    dbls[m++] = r[1];
}

//---------------------------------------------------------------------------//

void Vessel::track_LoR()
{
    Manager dmm(*pins[0]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);

    State_Vec_DMM states;
    states.initialize(1);

    vessel_kernel1<<<1,1>>>(
        pin, states.device_instance(), ints.data().get(), dbls.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    int n = 0, m = 0;
    double eps = 1.0e-10;

    EXPECT_TRUE(ints[n++]);

    EXPECT_EQ(State::VESSEL, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 1.2018173738e+01, eps);
    EXPECT_EQ(State::PLUS_X, ints[n++]);

    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 1.1616815398e+01, eps);
    EXPECT_EQ(State::MINUS_Y, ints[n++]);

    EXPECT_EQ(State::NONE, ints[n++]);
    EXPECT_EQ(State::VESSEL, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 1.4437568710e-02, eps);
    EXPECT_EQ(State::INTERNAL, ints[n++]);
    EXPECT_EQ(State::R0_VESSEL, ints[n++]);
    EXPECT_EQ(State::MODERATOR, ints[n++]);

    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 6.4107157490e+00, eps);
    EXPECT_EQ(State::INTERNAL, ints[n++]);
    EXPECT_EQ(State::R0_VESSEL, ints[n++]);
    EXPECT_EQ(State::VESSEL, ints[n++]);
}

//---------------------------------------------------------------------------//

void Vessel::track_HiR()
{
    Manager dmm(*pins[1]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);

    State_Vec_DMM states;
    states.initialize(1);

    vessel_kernel2<<<1,1>>>(
        pin, states.device_instance(), ints.data().get(), dbls.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    int n = 0, m = 0;
    double eps = 1.0e-10;

    EXPECT_TRUE(ints[n++]);

    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 1.2018173738e+01, eps);
    EXPECT_EQ(State::PLUS_X, ints[n++]);

    EXPECT_EQ(State::VESSEL, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 1.1616815398e+01, eps);
    EXPECT_EQ(State::MINUS_Y, ints[n++]);

    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_EQ(State::NONE, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 1.4437568710e-02, eps);
    EXPECT_EQ(State::INTERNAL, ints[n++]);
    EXPECT_EQ(State::R1_VESSEL, ints[n++]);
    EXPECT_EQ(State::VESSEL, ints[n++]);

    EXPECT_EQ(State::VESSEL, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 6.4107157490e+00, eps);
    EXPECT_EQ(State::INTERNAL, ints[n++]);
    EXPECT_EQ(State::R1_VESSEL, ints[n++]);
    EXPECT_EQ(State::MODERATOR, ints[n++]);
}

//---------------------------------------------------------------------------//

void Vessel::track_LoHiR()
{
    Manager dmm(*pins[2]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);

    State_Vec_DMM states;
    states.initialize(1);

    vessel_kernel3<<<1,1>>>(
        pin, states.device_instance(), ints.data().get(), dbls.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    int n = 0, m = 0;
    double eps = 1.0e-10;

    EXPECT_TRUE(ints[n++]);

    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 1.2018173738e+01, eps);
    EXPECT_EQ(State::PLUS_X, ints[n++]);

    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_EQ(State::NONE, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 1.4437568710e-02, eps);
    EXPECT_EQ(State::INTERNAL, ints[n++]);
    EXPECT_EQ(State::R1_VESSEL, ints[n++]);
    EXPECT_EQ(State::VESSEL, ints[n++]);

    EXPECT_EQ(State::VESSEL, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 3.4045225628e+00, eps);
    EXPECT_EQ(State::INTERNAL, ints[n++]);
    EXPECT_EQ(State::R0_VESSEL, ints[n++]);
    EXPECT_EQ(State::MODERATOR, ints[n++]);

    EXPECT_EQ(State::VESSEL, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 6.4107157490e+00, eps);
    EXPECT_EQ(State::INTERNAL, ints[n++]);
    EXPECT_EQ(State::R1_VESSEL, ints[n++]);
    EXPECT_EQ(State::MODERATOR, ints[n++]);
}

//---------------------------------------------------------------------------//

void Vessel::track_Hi2Lo()
{
    Manager dmm(*pins[2]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);

    State_Vec_DMM states;
    states.initialize(1);

    vessel_kernel4<<<1,1>>>(
        pin, states.device_instance(), ints.data().get(), dbls.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    int n = 0, m = 0;
    double eps = 1.0e-10;

    EXPECT_TRUE(ints[n++]);

    EXPECT_EQ(State::NONE, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++], 1.4437568710e-02, eps);
    EXPECT_EQ(State::INTERNAL, ints[n++]);
    EXPECT_EQ(State::R1_VESSEL, ints[n++]);
    EXPECT_EQ(State::VESSEL, ints[n++]);
    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_EQ(10, ints[n++]);

    EXPECT_EQ(State::R1_VESSEL, ints[n++]);
    EXPECT_EQ(State::VESSEL, ints[n++]);
    EXPECT_EQ(0, ints[n++]);
    EXPECT_EQ(State::INTERNAL, ints[n++]);

    EXPECT_SOFTEQ(dbls[m++], 5.3925058590e+00, eps);
    EXPECT_EQ(State::INTERNAL, ints[n++]);
    EXPECT_EQ(State::R0_VESSEL, ints[n++]);
    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_EQ(State::VESSEL, ints[n++]);
    EXPECT_EQ(101, ints[n++]);

    EXPECT_EQ(State::R0_VESSEL, ints[n++]);
    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_EQ(0, ints[n++]);
    EXPECT_EQ(State::INTERNAL, ints[n++]);

    EXPECT_SOFTEQ(dbls[m++], 1.0715916120e+01, eps);
    EXPECT_EQ(State::MINUS_X, ints[n++]);
    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_EQ(10, ints[n++]);

    EXPECT_SOFTEQ(-21.62*0.5, dbls[m++], eps);
}

//---------------------------------------------------------------------------//

void Vessel::track_Lo2Hi()
{
    Manager dmm(*pins[2]);

    auto pin = dmm.device_instance();

    thrust::device_vector<int>    ints(50, -1);
    thrust::device_vector<double> dbls(50, -1);

    State_Vec_DMM states;
    states.initialize(1);

    vessel_kernel5<<<1,1>>>(
        pin, states.device_instance(), ints.data().get(), dbls.data().get());

    thrust::host_vector<int>    rints(ints.begin(), ints.end());
    thrust::host_vector<double> rdbls(dbls.begin(), dbls.end());

    int n = 0, m = 0;
    double eps = 1.0e-10;

    EXPECT_TRUE(ints[n++]);

    EXPECT_EQ(State::NONE, ints[n++]);
    EXPECT_SOFTEQ(dbls[m++],  2.6917439842e+00, eps);
    EXPECT_EQ(State::INTERNAL, ints[n++]);
    EXPECT_EQ(State::R0_VESSEL, ints[n++]);
    EXPECT_EQ(State::VESSEL, ints[n++]);
    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_EQ(10, ints[n++]);

    EXPECT_EQ(State::R0_VESSEL, ints[n++]);
    EXPECT_EQ(State::VESSEL, ints[n++]);
    EXPECT_EQ(0, ints[n++]);
    EXPECT_EQ(State::INTERNAL, ints[n++]);

    EXPECT_SOFTEQ(dbls[m++], 5.1303903122e+00, eps);
    EXPECT_EQ(State::INTERNAL, ints[n++]);
    EXPECT_EQ(State::R1_VESSEL, ints[n++]);
    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_EQ(State::VESSEL, ints[n++]);
    EXPECT_EQ(101, ints[n++]);

    EXPECT_EQ(State::R1_VESSEL, ints[n++]);
    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_EQ(0, ints[n++]);
    EXPECT_EQ(State::INTERNAL, ints[n++]);

    EXPECT_SOFTEQ(dbls[m++], 1.4303925000e+01, eps);
    EXPECT_EQ(State::PLUS_Y, ints[n++]);
    EXPECT_EQ(State::MODERATOR, ints[n++]);
    EXPECT_EQ(10, ints[n++]);

    EXPECT_SOFTEQ(15.0, dbls[m++], eps);
}

//---------------------------------------------------------------------------//
// end of RTK_Cell_Tester_2.cu
//---------------------------------------------------------------------------//
