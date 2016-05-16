//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/Eigenvalue_Solver.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  Eigenvalue_Solver explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Eigenvalue_Solver.t.hh"
#include "solvers/LinAlgTypedefs.hh"

namespace profugus
{

template class Eigenvalue_Solver<EpetraTypes>;
template class Eigenvalue_Solver<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Eigenvalue_Solver.pt.cc
//---------------------------------------------------------------------------//
