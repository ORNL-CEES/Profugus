//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/test/Run_Launch_Args.hh
 * \author Steven Hamilton
 * \date   Tue Nov 24 14:09:29 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_test_Run_Launch_Args_hh
#define CudaUtils_test_Run_Launch_Args_hh

#include <vector>
#include "../cuda_utils/Host_Vector.hh"

//---------------------------------------------------------------------------//
// Launch_Args test functor.
template<typename Arch_T>
void run_launch_args(std::vector<double> &data);

#endif // CudaUtils_test_Run_Launch_Args_hh

//---------------------------------------------------------------------------//
//                        end of Run_Launch_Args.hh
//---------------------------------------------------------------------------//
