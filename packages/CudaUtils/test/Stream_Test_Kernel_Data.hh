//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Stream_Test_Kernel_Data.hh
 * \author Seth R Johnson
 * \date   Wed Oct 02 15:12:07 2013
 * \brief  Stream_Test_Kernel_Data kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Stream_Test_Kernel_Data_hh
#define cuda_utils_test_Stream_Test_Kernel_Data_hh

#include "../cuda_utils/Device_Vector.hh"
#include "../cuda_utils/Launch_Args.hh"

namespace cuda_utils
{
//---------------------------------------------------------------------------//

template<typename Arch_T, typename Float_T>
struct Stream_Test_Kernel_Data
{
    // >>> Typedefs
    typedef Device_Vector<Arch_T, Float_T> Device_Vector_Float;
    typedef Launch_Args<Arch_T>            Launch_Args_t;

    // >>> Data
    const int num_rays;
    const int num_segments;

    // Grid/block sizes
    Launch_Args_t       launch_args;

    // Stream/tau currently being swept
    Device_Vector_Float tau;

    Device_Vector_Float source;
    Device_Vector_Float input;
    Device_Vector_Float output;

    // Back-end buffers
    Device_Vector_Float tau_build;

    // >>> Constructor

    Stream_Test_Kernel_Data(
            unsigned int num_rays_,
            unsigned int num_segments_)
      : num_rays(num_rays_)
      , num_segments(num_segments_)
      , tau(num_rays * num_segments)
      , source(num_rays * num_segments)
      , input(num_rays)
      , output(num_rays)
      , tau_build(num_rays * num_segments)
    {
        /* * */
    }
};



//---------------------------------------------------------------------------//
} // end namespace cuda_utils
#endif // cuda_utils_test_Stream_Test_Kernel_Data_hh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Stream_Test_Kernel_Data.hh
//---------------------------------------------------------------------------//
