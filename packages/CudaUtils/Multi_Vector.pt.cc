//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Multi_Vector.pt.cc
 * \author Seth R Johnson
 * \date   Fri Aug 02 10:24:45 2013
 * \brief  Multi_Vector member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Multi_Vector.t.hh"

#include <cuda_utils/config.h>
#include "Device_Vector.hh"

namespace cuda
{
//---------------------------------------------------------------------------//

template class Multi_Vector<arch::Host, int>;
template class Multi_Vector<arch::Host, float>;
template class Multi_Vector<arch::Host, double>;

#ifdef USE_CUDA
template class Multi_Vector<arch::Device, int>;
template class Multi_Vector<arch::Device, float>;
template class Multi_Vector<arch::Device, double>;
#endif

//---------------------------------------------------------------------------//
} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of Multi_Vector.pt.cc
//---------------------------------------------------------------------------//
