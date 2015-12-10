//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/Run_Cuda_RNG.t.hh
 * \author Stuart Slattery
 * \date   Tue Nov 24 14:09:29 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Run_Cuda_RNG_t_hh
#define cuda_utils_test_Run_Cuda_RNG_t_hh

#include <memory>

#include "Run_Cuda_RNG.hh"

#include "../cuda_utils/Cuda_RNG.hh"
#include "../cuda_utils/Host_Vector.hh"

//---------------------------------------------------------------------------//
// Cuda_RNG test functor.
template<typename Arch>
void run_cuda_rng( cuda::Host_Vector<int>& seeds,
		   cuda::Host_Vector<cuda::Cuda_RNG>& rng,
		   cuda::Host_Vector<double>& data )
{
    // Copy input data to the device.
    std::shared_ptr<cuda::Device_Vector<Arch,int> > device_seeds =
	std::make_shared<cuda::Device_Vector<Arch,int> >( seeds );
    std::shared_ptr<cuda::Device_Vector<Arch,cuda::Cuda_RNG> > device_rng =
	std::make_shared<cuda::Device_Vector<Arch,cuda::Cuda_RNG> >( rng );
    std::shared_ptr<cuda::Device_Vector<Arch,double> > device_data =
	std::make_shared<cuda::Device_Vector<Arch,double> >( data );

    // Make launch args.
    cuda::Launch_Args<Arch> args;
    args.set_block_size(256);
    args.set_num_elements(data.size());

    // Create an rng init kernel.
    Cuda_RNG_Vector_Init_Kernel<Arch_T> 
	init_kernel( device_rng, device_seeds );

    // Create a vector fill kernel.
    Cuda_RNG_Vector_Fill_Kernel<Arch_T> 
	fill_kernel( device_rng, device_data );

    // Initialize the rngs.
    cuda::parallel_launch( init_kernel, args );

    // Fill the vector.
    cuda::parallel_launch( fill_kernel, args );

    // Copy the results back to the host.
    device_seeds->to_host( seeds );
    device_rng->to_host( rng );
    device_data->to_host( data );
}

#endif // cuda_utils_test_Run_Cuda_RNG_t_hh

//---------------------------------------------------------------------------//
//                        end of Run_Cuda_RNG.t.hh
//---------------------------------------------------------------------------//
