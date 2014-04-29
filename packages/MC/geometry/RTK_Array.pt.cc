//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/RTK_Array.pt.cc
 * \author Thomas M. Evans
 * \date   Tue Dec 14 12:38:39 2010
 * \brief  RTK_Array explicit instantiations.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RTK_Cell.hh"
#include "RTK_Array.t.hh"

namespace denovo
{

template class RTK_Array<RTK_Cell>;
template class RTK_Array< RTK_Array<RTK_Cell> >;
template class RTK_Array< RTK_Array< RTK_Array<RTK_Cell> > >;

} // end namespace denovo

//---------------------------------------------------------------------------//
//                 end of RTK_Array.pt.cc
//---------------------------------------------------------------------------//
