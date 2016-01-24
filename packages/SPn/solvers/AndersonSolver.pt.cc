//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/AndersonSolver.pt.cc
 * \author Steven Hamilton
 * \date   Wed Apr 01 11:01:28 2015
 * \brief  AndersonSolver explicit instantiations.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "AndersonSolver.t.hh"

#include "LinAlgTypedefs.hh"

namespace profugus
{

template class AndersonSolver<EpetraTypes>;
template class AndersonSolver<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of AndersonSolver.pt.cc
//---------------------------------------------------------------------------//
