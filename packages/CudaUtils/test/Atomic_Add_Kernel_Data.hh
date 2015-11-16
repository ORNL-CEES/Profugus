//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/test/Atomic_Add_Kernel_Data.hh
 * \author Seth R Johnson
 * \date   Thu Aug 15 11:11:28 2013
 * \brief  Atomic_Add_Kernel_Data kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Atomic_Add_Kernel_Data_hh
#define cuda_utils_test_Atomic_Add_Kernel_Data_hh

#include "../Vector_Traits.hh"
#include "../Device_Vector.hh"
#include "../Launch_Args.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
//! Lock kernel data wrapper
template<typename Arch_Switch, typename Float_T>
struct Atomic_Add_Kernel_Data
{
    // >>> Typedefs
    typedef cuda::Vector_Traits<Arch_Switch, Float_T> Traits_t;
    typedef typename Traits_t::Device_Vector_Float    Device_Vector_Float;
    typedef Launch_Args<Arch_Switch>                  Launch_Args_t;

    // >>> Data
    Launch_Args_t launch_args;
    unsigned int  num_increments;

    Device_Vector_Float output;

    // >>> Constructor

    Atomic_Add_Kernel_Data()
      : num_increments(0)
      , output(1)
    {
        /* * */
    }
};

//---------------------------------------------------------------------------//
} // end namespace cuda

#endif // cuda_utils_test_Atomic_Add_Kernel_Data_hh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/test/Atomic_Add_Kernel_Data.hh
//---------------------------------------------------------------------------//
