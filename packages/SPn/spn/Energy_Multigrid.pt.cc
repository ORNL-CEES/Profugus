//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Energy_Multigrid.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  Energy_Multigrid explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Energy_Multigrid.t.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"

namespace profugus
{

template class Energy_Multigrid<EpetraTypes>;
template class Energy_Multigrid<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Energy_Multigrid.pt.cc
//---------------------------------------------------------------------------//
