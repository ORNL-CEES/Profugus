//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/Launch_Args_Kernel.pt.cc
 * \author Stuart Slattery
 * \date   Tue Nov 24 14:09:29 2015
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../cuda_utils/Definitions.hh"
#include "../cuda_utils/Pseudo_Cuda.hh"
#include "../cuda_utils/Launch_Args.t.hh"

#include "Launch_Args_Kernel.hh"

typedef cuda::arch::Host Host;
template class Functor<Host>;
template void cuda::parallel_launch<Functor<Host> >(
        Functor<Host> &, const cuda::Launch_Args<Host> &);

//---------------------------------------------------------------------------//
//                        end of Launch_Args_Kernel.pt.cc
//---------------------------------------------------------------------------//
