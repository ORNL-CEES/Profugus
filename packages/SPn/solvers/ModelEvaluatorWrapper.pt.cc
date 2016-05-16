//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/ModelEvaluatorWrapper.pt.cc
 * \author Steven Hamilton
 * \date   Wed Apr 01 12:39:40 2015
 * \brief  ModelEvaluatorWrapper explicit instantiations
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "ModelEvaluatorWrapper.t.hh"

#include "LinAlgTypedefs.hh"

namespace profugus
{

template class ModelEvaluatorWrapper<EpetraTypes>;
template class ModelEvaluatorWrapper<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of ModelEvaluatorWrapper.pt.cc
//---------------------------------------------------------------------------//
