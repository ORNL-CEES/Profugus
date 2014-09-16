//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/PreconditionerBuilder.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  PreconditionerBuilder explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "PreconditionerBuilder.t.hh"

#include "LinAlgTypedefs.hh"

namespace profugus
{

template class PreconditionerBuilder<EpetraTypes>;
template class PreconditionerBuilder<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of PreconditionerBuilder.pt.cc
//---------------------------------------------------------------------------//
