//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Linear_System.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  Linear_System explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Linear_System.t.hh"
#include "solvers/LinAlgTypedefs.hh"

namespace profugus
{

template class Linear_System<EpetraTypes>;
//template class Linear_System<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Linear_System.pt.cc
//---------------------------------------------------------------------------//
