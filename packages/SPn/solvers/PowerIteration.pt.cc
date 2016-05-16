//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/PowerIteration.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  PowerIteration explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "PowerIteration.t.hh"

#include "LinAlgTypedefs.hh"

namespace profugus
{

template class PowerIteration<EpetraTypes>;
template class PowerIteration<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of PowerIteration.pt.cc
//---------------------------------------------------------------------------//
