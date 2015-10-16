//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Particle_Metaclass.pt.cc
 * \author Thomas M. Evans
 * \date   Wed Jul 23 23:46:36 2014
 * \brief  Instantiate metaclass on Particle.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "utils/Metaclass.t.hh"
#include "Particle.hh"

namespace profugus
{

template class Metaclass<Particle>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Particle_Metaclass.pt.cc
//---------------------------------------------------------------------------//
