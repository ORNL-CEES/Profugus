//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/StratimikosSolver.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  StratimikosSolver explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "StratimikosSolver.t.hh"

#include "TpetraTypedefs.hh"

namespace profugus
{

template class StratimikosSolver<Epetra_MultiVector,Epetra_Operator>;
template class StratimikosSolver<Tpetra_MultiVector,Tpetra_Operator>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of StratimikosSolver.pt.cc
//---------------------------------------------------------------------------//
