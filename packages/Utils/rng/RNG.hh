//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/RNG.hh
 * \author Thomas M. Evans
 * \date   Fri Apr 25 14:57:29 2014
 * \brief  RNG class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef rng_RNG_hh
#define rng_RNG_hh

#include <vector>
#include <random>

#include <Utils/config.h>
#include "harness/DBC.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class RNG
 * \brief Random-number generator class interface.
 */
/*!
 * \example rng/test/tstRNG.cc
 *
 * Test of RNG.
 */
//===========================================================================//

class RNG
{
  private:
    // >>> DATA

    // random engine.
    std::mt19937_64 d_engine;

    // distribution.
    std::uniform_real_distribution<double> d_dist;

  public:
    // Constructors
    inline RNG() : d_engine(0) {}
    inline RNG(int seed) : d_engine(seed) {}

    // >>> Services provided by RNG class.

    //! Get a uniform random number on [0,1] (closed on both sides)
    double ran() { return d_dist(d_engine); }
};

} // end namespace profugus

#endif // rng_RNG_hh

//---------------------------------------------------------------------------//
//                 end of RNG.hh
//---------------------------------------------------------------------------//
