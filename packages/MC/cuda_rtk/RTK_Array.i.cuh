//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_Array.i.cuh
 * \author Tom Evans
 * \date   Wed Jan 04 15:43:43 2017
 * \brief  RTK_Array CUDA device class definitions.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_RTK_Array_i_cuh
#define MC_cuda_rtk_RTK_Array_i_cuh

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Find the object a point is in.
 */
template<class T>
__device__
int RTK_Array<T>::find_object(
    const Space_Vector &r,
    Geo_State_t        &state) const
{
    return 0;
}

} // end namespace cuda_profugus

#endif // MC_cuda_rtk_RTK_Array_i_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Array.i.cuh
//---------------------------------------------------------------------------//
