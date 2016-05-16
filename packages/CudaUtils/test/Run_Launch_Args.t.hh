//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/test/Run_Launch_Args.t.hh
 * \author Steven Hamilton
 * \date   Tue Nov 24 14:09:29 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_test_Run_Launch_Args_t_hh
#define CudaUtils_test_Run_Launch_Args_t_hh

#include "Run_Launch_Args.hh"

#include "Launch_Args_Kernel.hh"
#include "../cuda_utils/Launch_Args.hh"
#include "../cuda_utils/Host_Vector.hh"

//---------------------------------------------------------------------------//
// Launch_Args test functor.
template<typename Arch>
void run_launch_args(std::vector<double> &data)
{
    // Make launch args.
    cuda::Launch_Args<Arch> args;
    args.set_block_size(256);
    args.set_num_elements(1024);

    // Create a functor.
    double value = 2.3;
    int size = args.num_elements();
    Functor<Arch> functor( size, value );

    // Call the kernel launch.
    cuda::parallel_launch( functor, args );

    // Get the data from the functor.
    cuda::Host_Vector<double> host_vec(size);
    functor.assign_data(host_vec);

    // Copy data into function argument
    data.clear();
    data.insert( data.end(), host_vec.begin(), host_vec.end() );
}

#endif // CudaUtils_test_Run_Launch_Args_t_hh

//---------------------------------------------------------------------------//
//                        end of Run_Launch_Args.t.hh
//---------------------------------------------------------------------------//
