//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/RayleighQuotient.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  RayleighQuotient explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "RayleighQuotient.t.hh"

#include "TpetraTypedefs.hh"

namespace profugus
{

template class RayleighQuotient<Epetra_MultiVector,Epetra_Operator>;
template class RayleighQuotient<Tpetra_MultiVector,Tpetra_Operator>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of RayleighQuotient.pt.cc
//---------------------------------------------------------------------------//
