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

namespace profugus
{

template class RayleighQuotient<Epetra_MultiVector,Epetra_Operator>;

typedef KokkosClassic::SerialNode Node;
typedef Tpetra::MultiVector<double,int,int,Node> MV;
typedef Tpetra::Operator<double,int,int,Node> OP;
template class RayleighQuotient<MV,OP>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of RayleighQuotient.pt.cc
//---------------------------------------------------------------------------//
