//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Arnoldi.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  Arnoldi explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "Arnoldi.t.hh"

#include "TpetraTypedefs.hh"

namespace profugus
{

template class Arnoldi<Epetra_MultiVector,Epetra_Operator>;
template class Arnoldi<Tpetra_MultiVector,Tpetra_Operator>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Arnoldi.pt.cc
//---------------------------------------------------------------------------//
