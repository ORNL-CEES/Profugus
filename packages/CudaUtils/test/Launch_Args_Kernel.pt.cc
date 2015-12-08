//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/Launch_Args_Kernel.cc
 * \author Stuart Slattery
 * \date   Tue Nov 24 14:09:29 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../cuda_utils/Definitions.hh"
#include "../cuda_utils/Pseudo_Cuda.hh"

#include "Launch_Args_Kernel.hh"

template<>
void Functor<cuda::arch::Host>::operator()(const std::size_t idx )
{
    d_device_data[idx] += idx;
}

//---------------------------------------------------------------------------//
//                        end of Launch_Args_Kernel.cc
//---------------------------------------------------------------------------//
