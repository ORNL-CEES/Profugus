//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/Solver_Base.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  Solver_Base explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Solver_Base.hh"
#include "solvers/LinAlgTypedefs.hh"

namespace profugus
{


template class Solver_Base_Tmpl<EpetraTypes>;
template class Solver_Base_Tmpl<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Solver_Base.pt.cc
//---------------------------------------------------------------------------//
