//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/Polyglot_Kernel_Data.hh
 * \author Seth R Johnson
 * \date   Wed Aug 14 13:03:02 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_test_Polyglot_Kernel_Data_hh
#define cuda_utils_test_Polyglot_Kernel_Data_hh

#include "../Vector_Traits.hh"
#include "../Device_Vector.hh"
#include "../Launch_Args.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
//! Polyglot kernel data wrapper
template<typename Arch_Switch>
struct Polyglot_Kernel_Data
{
    // >>> Typedefs
    typedef cuda::Vector_Traits<Arch_Switch, float> Traits_t;
    typedef typename Traits_t::Device_Vector_Float  Device_Vector_Float;
    typedef Launch_Args<Arch_Switch>                Launch_Args_t;

    // >>> Data
    Launch_Args_t launch_args;

    Device_Vector_Float input;
    Device_Vector_Float output;

    // >>> Constructor

    Polyglot_Kernel_Data(unsigned int num_elements)
      : input(num_elements)
      , output(num_elements)
    {
        /* * */
    }
};

//---------------------------------------------------------------------------//
} // end namespace cuda

#endif // cuda_utils_test_Polyglot_Kernel_Data_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Polyglot_Kernel_Data.hh
//---------------------------------------------------------------------------//
