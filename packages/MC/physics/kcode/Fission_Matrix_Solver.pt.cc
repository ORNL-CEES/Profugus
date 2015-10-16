//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Matrix_Solver.pt.cc
 * \author Thomas M. Evans
 * \date   Mon Dec 08 17:18:34 2014
 * \brief  Fission_Matrix_Solver member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "solvers/LinAlgTypedefs.hh"
#include "Fission_Matrix_Solver.t.hh"

namespace profugus
{

template class Fission_Matrix_Solver<EpetraTypes>;
template class Fission_Matrix_Solver<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Solver.pt.cc
//---------------------------------------------------------------------------//
