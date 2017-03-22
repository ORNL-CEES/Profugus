//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_State_Tester.cu
 * \author Tom Evans
 * \date   Fri Nov 18 15:25:34 2016
 * \brief  RTK_State_Tester kernel definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "Utils/gtest/Gtest_Functions.hh"
#include "../RTK_State.cuh"

//---------------------------------------------------------------------------//
// TYPES
//---------------------------------------------------------------------------//

using State        = cuda_profugus::RTK_State;
using Space_Vector = State::Space_Vector;
using Coordinates  = State::Coordinates;

//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//

__global__
void on_device(Space_Vector *r,
               Space_Vector *dir,
               Coordinates  *c,
               int          *i)
{
    State state;

    state.d_r   = {0.1, 0.2, 0.3};
    state.d_dir = {1.0/std::sqrt(3.0),
                   1.0/std::sqrt(3.0),
                   1.0/std::sqrt(3.0)};

    state.level_coord[0] = {1, 1, 1};
    state.level_coord[1] = {5, 4, 0};
    state.level_coord[2] = {12, 3, 7};

    *r   = state.d_r;
    *dir = state.d_dir;

    c[0] = state.level_coord[0];
    c[1] = state.level_coord[1];
    c[2] = state.level_coord[2];

    state.escaping_face = State::MINUS_Z;
    state.region        = 2;
    state.face          = State::PLUS_X;

    i[0] = state.escaping_face;
    i[1] = state.region;
    i[2] = state.face;
    i[3] = State::plus_face(1);
}

//---------------------------------------------------------------------------//

__global__
void to_device(State *s)
{
    State &s0 = s[0];
    State &s1 = s[1];
    State &s2 = s[2];

    s0.d_r = {0.1, 0.2, 0.3};
    s1.d_r = {21.1, 5.3, 4.1};
    s2.d_r = {3.1, 6.0, 1.0};
}

//---------------------------------------------------------------------------//
// HOST TEST FUNCTIONS
//---------------------------------------------------------------------------//

namespace rtk_state_test
{

//---------------------------------------------------------------------------//

void test_on_device()
{
    thrust::device_vector<Space_Vector> r(1);
    thrust::device_vector<Space_Vector> d(1);
    thrust::device_vector<Coordinates>  c(3);
    thrust::device_vector<int>          i(4);

    on_device<<<1,1>>>(r.data().get(),
                       d.data().get(),
                       c.data().get(),
                       i.data().get());

    thrust::host_vector<Space_Vector> ret_r(r.begin(), r.end());
    thrust::host_vector<Space_Vector> ret_d(d.begin(), d.end());
    thrust::host_vector<Coordinates>  ret_c(c.begin(), c.end());
    thrust::host_vector<int>          ret_i(i.begin(), i.end());


    EXPECT_EQ(0.1, ret_r[0][0]);
    EXPECT_EQ(0.2, ret_r[0][1]);
    EXPECT_EQ(0.3, ret_r[0][2]);

    EXPECT_SOFT_EQ(1.0/std::sqrt(3.0), ret_d[0][0]);
    EXPECT_SOFT_EQ(1.0/std::sqrt(3.0), ret_d[0][1]);
    EXPECT_SOFT_EQ(1.0/std::sqrt(3.0), ret_d[0][2]);

    EXPECT_EQ(1, ret_c[0][0]);
    EXPECT_EQ(1, ret_c[0][1]);
    EXPECT_EQ(1, ret_c[0][2]);

    EXPECT_EQ(5, ret_c[1][0]);
    EXPECT_EQ(4, ret_c[1][1]);
    EXPECT_EQ(0, ret_c[1][2]);

    EXPECT_EQ(12, ret_c[2][0]);
    EXPECT_EQ(3,  ret_c[2][1]);
    EXPECT_EQ(7,  ret_c[2][2]);

    EXPECT_EQ(1004, ret_i[0]);
    EXPECT_EQ(2,    ret_i[1]);
    EXPECT_EQ(1001, ret_i[2]);
    EXPECT_EQ(State::PLUS_Y, ret_i[3]);
}

//---------------------------------------------------------------------------//

void test_to_device()
{
    thrust::device_vector<State> s(3);

    to_device<<<1,1>>>(s.data().get());

    thrust::host_vector<State> result(s.begin(), s.end());

    EXPECT_EQ(0.1, result[0].d_r[0]);
    EXPECT_EQ(0.2, result[0].d_r[1]);
    EXPECT_EQ(0.3, result[0].d_r[2]);

    EXPECT_EQ(21.1, result[1].d_r[0]);
    EXPECT_EQ(5.3,  result[1].d_r[1]);
    EXPECT_EQ(4.1,  result[1].d_r[2]);

    EXPECT_EQ(3.1, result[2].d_r[0]);
    EXPECT_EQ(6.0, result[2].d_r[1]);
    EXPECT_EQ(1.0, result[2].d_r[2]);
}

} // end namespace rtk_state_test

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/RTK_State_Tester.cu
//---------------------------------------------------------------------------//
