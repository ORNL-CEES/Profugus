//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Linear_System_FV.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  Linear_System_FV explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Linear_System_FV.t.hh"
#include "solvers/LinAlgTypedefs.hh"

namespace profugus
{

template class Linear_System_FV<EpetraTypes>;
//template class Linear_System_FV<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Linear_System_FV.pt.cc
//---------------------------------------------------------------------------//
