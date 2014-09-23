//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/InverseOperator.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  InverseOperator explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "InverseOperator.hh"

#include "LinAlgTypedefs.hh"

namespace profugus
{

template class InverseOperator<EpetraTypes>;
template class InverseOperator<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of InverseOperator.pt.cc
//---------------------------------------------------------------------------//
