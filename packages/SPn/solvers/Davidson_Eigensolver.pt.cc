//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/Davidson_Eigensolver.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  Davidson_Eigensolver explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "Davidson_Eigensolver.t.hh"

#include "LinAlgTypedefs.hh"

namespace profugus
{

template class Davidson_Eigensolver<EpetraTypes>;
template class Davidson_Eigensolver<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Davidson_Eigensolver.pt.cc
//---------------------------------------------------------------------------//
