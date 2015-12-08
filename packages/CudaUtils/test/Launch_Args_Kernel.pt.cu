//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/Launch_Args_Kernel.cu
 * \author Stuart Slattery
 * \date   Tue Nov 24 14:09:29 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Launch_Args_Kernel.hh"

#include "../cuda_utils/Definitions.hh"

__device__
template<>
void Functor<cuda::arch::Device>::operator()(const std::size_t idx )
{
    d_device_data[idx] += idx;
}

//---------------------------------------------------------------------------//
//                        end of Launch_Args_Kernel.cu
//---------------------------------------------------------------------------//
