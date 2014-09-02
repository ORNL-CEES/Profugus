//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/ShiftedOperator.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  ShiftedOperator explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "ShiftedOperator.hh"

#include "TpetraTypedefs.hh"

namespace profugus
{

template class ShiftedOperator<Epetra_MultiVector,Epetra_Operator>;
template class ShiftedOperator<Tpetra_MultiVector,Tpetra_Operator>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of ShiftedOperator.pt.cc
//---------------------------------------------------------------------------//
