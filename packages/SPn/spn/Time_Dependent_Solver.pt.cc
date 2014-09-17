//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Time_Dependent_Solver.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  Time_Dependent_Solver explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Time_Dependent_Solver.t.hh"
#include "solvers/LinAlgTypedefs.hh"

namespace profugus
{

template class Time_Dependent_Solver<EpetraTypes>;
//template class Time_Dependent_Solver<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Time_Dependent_Solver.pt.cc
//---------------------------------------------------------------------------//
