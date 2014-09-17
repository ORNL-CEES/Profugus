//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Fixed_Source_Solver.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  Fixed_Source_Solver explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fixed_Source_Solver.t.hh"
#include "solvers/LinAlgTypedefs.hh"

namespace profugus
{

template class Fixed_Source_Solver<EpetraTypes>;
//template class Fixed_Source_Solver<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.pt.cc
//---------------------------------------------------------------------------//
