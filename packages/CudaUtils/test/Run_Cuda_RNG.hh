//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/Run_Cuda_RNG.hh
 * \author Stuart Slattery
 * \date   Tue Nov 24 14:09:29 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Run_Cuda_RNG_hh
#define cuda_utils_test_Run_Cuda_RNG_hh

#include <vector>
#include "../cuda_utils/Host_Vector.hh"
#include "../cuda_utils/Cuda_RNG.hh"

//---------------------------------------------------------------------------//
// Cuda_RNG test function.
template<typename Arch_T>
void run_cuda_rng(cuda::Host_Vector<int>& seeds,
		  cuda::Host_Vector<cuda::Cuda_RNG>& rng,
		  cuda::Host_Vector<double> &data);

#endif // cuda_utils_test_Run_Cuda_RNG_hh

//---------------------------------------------------------------------------//
//                        end of Run_Cuda_RNG.hh
//---------------------------------------------------------------------------//
