//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Anderson_Solver.pt.cc
 * \author Thomas M. Evans
 * \date   Tue Apr 07 20:53:11 2015
 * \brief  Anderson_Solver explicit instantiations.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "solvers/LinAlgTypedefs.hh"
#include "Anderson_Solver.t.hh"

namespace profugus
{

template class Anderson_Solver<EpetraTypes>;
template class Anderson_Solver<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
// end of MC/mc/Anderson_Solver.pt.cc
//---------------------------------------------------------------------------//
