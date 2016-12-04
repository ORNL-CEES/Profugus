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

using Device_Cell  = cuda_profugus::RTK_Cell;
using DMM          = cuda_profugus::RTK_Cell_DMM;
using Space_Vector = Device_Cell::Space_Vector;
using State        = Device_Cell::Geo_State_t;

using def::X; using def::Y; using def::Z;

//---------------------------------------------------------------------------//
// SINGLE-SHELL
//---------------------------------------------------------------------------//

__global__
void single_shell_kernel1(
    Device_Cell   pin,
    int          *ints,
    double	 *dbls,
    Space_Vector *svs)
{
    auto r  = pin.radii();
    dbls[0] = r[0];
    dbls[1] = pin.pitch(X);
    dbls[2] = pin.pitch(Y);

    ints[0] = pin.num_regions();
    ints[1] = pin.num_shells();
    //ints[2] = pin.region(0.0, 0.5401);
    //ints[3] = pin.region(0.0, 0.5399);
    ints[4] = pin.matid(0);
    ints[5] = pin.matid(1);

    pin.get_extents(svs[0], svs[1]);
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

    thrust::device_vector<int>          ints(20, -1);
    thrust::device_vector<double>       dbls(20, -1);
    thrust::device_vector<Space_Vector> svs(20);

    {

        single_shell_kernel1<<<1,1>>>(p1, ints.data().get(),
dbls.data().get(), svs.data().get() );

        thrust::host_vector<int>          rints(ints.begin(), ints.end());
        thrust::host_vector<double>       rdbls(dbls.begin(), dbls.end());
        thrust::host_vector<Space_Vector> rsvs(svs.begin(), svs.end());

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
    }

    {
        single_shell_kernel2<<<1,1>>>(p2);
    }
}

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/RTK_Cell_Tester.cu
//---------------------------------------------------------------------------//
