//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Energy_Prolongation.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  Energy_Prolongation explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Energy_Prolongation.t.hh"
#include "solvers/LinAlgTypedefs.hh"

namespace profugus
{

template class Energy_Prolongation<EpetraTypes>;
//template class Energy_Prolongation<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Energy_Prolongation.pt.cc
//---------------------------------------------------------------------------//
