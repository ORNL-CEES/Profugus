//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/RTK_Geometry.pt.cc
 * \author Thomas M. Evans
 * \date   Tue Jan 25 10:02:33 2011
 * \brief  Explicit instantiations of RTK_Geometry types.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Geometry.t.hh"

namespace profugus
{

template class Geometry< RTK_Array<RTK_Cell> >;
template class Geometry< RTK_Array< RTK_Array<RTK_Cell> > >;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of RTK_Geometry.pt.cc
//---------------------------------------------------------------------------//
