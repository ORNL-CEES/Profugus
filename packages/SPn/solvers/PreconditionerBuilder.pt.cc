//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/PreconditionerBuilder.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  PreconditionerBuilder explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Tpetra_Operator.hpp"
#include "Tpetra_MultiVector.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "PreconditionerBuilder.t.hh"

namespace profugus
{

template class PreconditionerBuilder<Epetra_Operator>;

typedef KokkosClassic::SerialNode Node;
typedef Tpetra::Operator<double,int,int,Node> OP;
template class PreconditionerBuilder<OP>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of PreconditionerBuilder.pt.cc
//---------------------------------------------------------------------------//
