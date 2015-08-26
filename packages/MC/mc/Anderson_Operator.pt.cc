//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Anderson_Operator.pt.cc
 * \author Thomas M. Evans
 * \date   Tue Apr 07 21:05:56 2015
 * \brief  Anderson_Operator explicit instantiations.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "solvers/LinAlgTypedefs.hh"
#include "Anderson_Operator.t.hh"

namespace profugus
{

template class Anderson_Operator<EpetraTypes>;
template class Anderson_Operator<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
// end of MC/mc/Anderson_Operator.pt.cc
//---------------------------------------------------------------------------//
