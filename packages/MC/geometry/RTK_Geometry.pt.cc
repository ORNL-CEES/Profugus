//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/rtk/RTK_Geometry.pt.cc
 * \author Thomas M. Evans
 * \date   Tue Jan 25 10:02:33 2011
 * \brief  Explicit instantiations of RTK_Geometry types.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RTK_Geometry.t.hh"

namespace denovo
{

template class RTK_Geometry< RTK_Array<RTK_Cell> >;
template class RTK_Geometry< RTK_Array< RTK_Array<RTK_Cell> > >;

} // end namespace denovo

//---------------------------------------------------------------------------//
//                 end of RTK_Geometry.pt.cc
//---------------------------------------------------------------------------//
