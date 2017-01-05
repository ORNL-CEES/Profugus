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
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class RTK_Array<RTK_Cell>;
template class RTK_Array< RTK_Array<RTK_Cell> >;

//---------------------------------------------------------------------------//
// RTK_ARRAY_DMM MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
RTK_Core_Array_DMM::RTK_Core_Array_DMM(
    const Host_Core_Array &host_array)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a device instance.
 */
RTK_Core_Array_DMM::Core_Array
RTK_Core_Array_DMM::device_instance()
{
    return Core_Array();
}

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Array.cu
//---------------------------------------------------------------------------//
