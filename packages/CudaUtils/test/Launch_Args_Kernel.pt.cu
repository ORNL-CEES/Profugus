//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/Launch_Args_Kernel.pt.cu
 * \author Stuart Slattery
 * \date   Tue Nov 24 14:09:29 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Launch_Args_Kernel.hh"

#include "../cuda_utils/Definitions.hh"
#include "../cuda_utils/Launch_Args.t.cuh"

typedef cuda::arch::Device Device;
template class Functor<Device>;
template void cuda::parallel_launch<Functor<Device> >(
        Functor<Device> &, const cuda::Launch_Args<Device> &);

//---------------------------------------------------------------------------//
//                        end of Launch_Args_Kernel.pt.cu
//---------------------------------------------------------------------------//
