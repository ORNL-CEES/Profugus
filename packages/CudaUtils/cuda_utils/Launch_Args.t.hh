//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Launch_Args.t.hh
 * \author Stuart Slattery
 * \date   Wed Oct 02 13:16:37 2013
 * \brief  Launch class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Launch_Args_t_hh
#define cuda_utils_Launch_Args_t_hh

#include "Pseudo_Cuda.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// Host specialization.
template<class Kernel>
void ParallelLaunch<cuda::arch::Host>::launch(
    Kernel& kernel, const Launch_Args<cuda::arch::Host>& launch_args )
{
    REQUIRE( launch_args.is_valid() );
    std::size_t num_t = launch_args.grid_size*launch_args.block_size;
    for ( std::size_t n = 0; n < num_t; ++n )
    {
	kernel( n );
    }
}

//---------------------------------------------------------------------------//

} // end namespace cuda

#endif // cuda_utils_Launch_Args_t_hh

//---------------------------------------------------------------------------//
//                 end of Launch_Args.t.hh
//---------------------------------------------------------------------------//
