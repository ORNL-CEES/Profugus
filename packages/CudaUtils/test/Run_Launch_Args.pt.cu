//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/Launch_Args_Kernel.pt.cu
 * \author Stuart Slattery
 * \date   Tue Nov 24 14:09:29 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../cuda_utils/Definitions.hh"
#include "../cuda_utils/Launch_Args.t.cuh"
#include "Run_Launch_Args.t.hh"

typedef cuda_utils::arch::Device Device;
template void run_launch_args<Device>(std::vector<double> &);

//---------------------------------------------------------------------------//
//                        end of Launch_Args_Kernel.pt.cu
//---------------------------------------------------------------------------//
