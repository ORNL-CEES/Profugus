//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Solver.t.cuh
 * \author Thomas M. Evans
 * \date   Tue May 13 14:56:00 2014
 * \brief  Solver member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Solver_t_hh
#define mc_Solver_t_hh

#include "Solver.hh"

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Virtual destructor.
 */
template <class Geometry>
Solver<Geometry>::~Solver()
{
}

} // end namespace cuda_profugus

#endif // mc_Solver_t_hh

//---------------------------------------------------------------------------//
//                 end of Solver.t.hh
//---------------------------------------------------------------------------//
