//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Host_Vector.pt.cc
 * \author Seth R Johnson
 * \date   Mon Aug 12 08:48:53 2013
 * \brief  Host_Vector template instantiations
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Host_Vector.t.hh"
#include "Definitions.hh"

namespace cuda
{
typedef cuda_utils::Space_Vector Space_Vector;

template class Host_Vector<int>;
template class Host_Vector<unsigned int>;
template class Host_Vector<float>;
template class Host_Vector<double>;

} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of Host_Vector.pt.cc
//---------------------------------------------------------------------------//
