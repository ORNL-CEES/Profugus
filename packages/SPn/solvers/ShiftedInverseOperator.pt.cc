//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/ShiftedInverseOperator.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  ShiftedInverseOperator explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "ShiftedInverseOperator.hh"

#include "LinAlgTypedefs.hh"

namespace profugus
{

template class ShiftedInverseOperator<EpetraTypes>;
template class ShiftedInverseOperator<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of ShiftedInverseOperator.pt.cc
//---------------------------------------------------------------------------//
