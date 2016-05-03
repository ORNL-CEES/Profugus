//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Device_Vector_Lite.hh
 * \author Steven Hamilton
 * \date   Mon May 02 15:19:21 2016
 * \brief  Device_Vector_Lite class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Device_Vector_Lite_hh
#define CudaUtils_cuda_utils_Device_Vector_Lite_hh

#include "CudaDBC.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Device_Vector_Lite
 * \brief On-device implementation of Vector_Lite
 */
/*!
 * \example cuda_utils/test/tstDevice_Vector_Lite.cc
 *
 * Test of Device_Vector_Lite.
 */
//===========================================================================//

template <class T, size_t N>
class Device_Vector_Lite
{
  public:

    __host__ __device__ T & operator[](int i)
    {
        REQUIRE( i < N );
        return d_data[i];
    }

    __host__ __device__ const T & operator[](int i) const
    {
        REQUIRE( i < N );
        return d_data[i];
    }

    // This data member is public so that the class satisfies the
    // C++11 is_standard_layout concept.
    T d_data[N];
};

//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
#endif // CudaUtils_cuda_utils_Device_Vector_Lite_hh

//---------------------------------------------------------------------------//
// end of CudaUtils/cuda_utils/Device_Vector_Lite.hh
//---------------------------------------------------------------------------//
