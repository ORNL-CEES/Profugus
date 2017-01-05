//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_Array.cu
 * \author Tom Evans
 * \date   Wed Jan 04 15:43:43 2017
 * \brief  RTK_Array member and kernel definitions.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RTK_Array.cuh"
#include "RTK_Cell.cuh"

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
// RTK_ARRAY MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class T>
RTK_Array<T>::RTK_Array()
{
}

//---------------------------------------------------------------------------//
// RTK_ARRAY_DMM MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class T>
RTK_Array_DMM<T>::RTK_Array_DMM(
    const Host_RTK_Array &host_array)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a device instance.
 */
template<class T>
RTK_Array<T> RTK_Array_DMM<T>::device_instance()
{
    return RTK_Array<T>();
}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class RTK_Array<RTK_Cell>;
template class RTK_Array< RTK_Array<RTK_Cell> >;

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Array.cu
//---------------------------------------------------------------------------//
