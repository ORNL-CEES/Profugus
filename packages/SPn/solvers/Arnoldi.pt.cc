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

#include "LinAlgTypedefs.hh"

namespace profugus
{

template class Arnoldi<EpetraTypes>;
template class Arnoldi<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Arnoldi.pt.cc
//---------------------------------------------------------------------------//
