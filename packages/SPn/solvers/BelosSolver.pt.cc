//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/BelosSolver.pt.cc
 * \author Steven Hamilton
 * \brief  BelosSolver explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "BelosEpetraAdapter.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosSolver.t.hh"

#include "LinAlgTypedefs.hh"

namespace profugus
{

template class BelosSolver<EpetraTypes>;
template class BelosSolver<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of BelosSolver.pt.cc
//---------------------------------------------------------------------------//
