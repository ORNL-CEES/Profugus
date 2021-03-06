//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Launch_Args.t.hh
 * \author Stuart Slattery
 * \date   Wed Oct 02 13:16:37 2013
 * \brief  Launch class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Launch_Args_t_hh
#define CudaUtils_cuda_utils_Launch_Args_t_hh

#include "harness/DBC.hh"
#include "Launch_Args.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// Host specialization.
template <class Kernel>
void parallel_launch(
    Kernel& kernel, const Launch_Args<cuda::arch::Host>& launch_args )
{
    REQUIRE( launch_args.is_valid() );
    std::size_t num_t = launch_args.num_elements();
    for ( std::size_t n = 0; n < num_t; ++n )
    {
        kernel( n );
    }
}

//---------------------------------------------------------------------------//

} // end namespace cuda

#endif // CudaUtils_cuda_utils_Launch_Args_t_hh

//---------------------------------------------------------------------------//
//                 end of Launch_Args.t.hh
//---------------------------------------------------------------------------//
