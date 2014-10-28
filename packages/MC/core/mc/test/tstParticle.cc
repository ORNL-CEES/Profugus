//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstParticle.cc
 * \author Thomas M. Evans
 * \date   Fri Apr 25 11:32:36 2014
 * \brief  Particle class test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Particle.hh"

#include "gtest/utils_gtest.hh"

using profugus::Particle;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Particle, construction)
{
    Particle p;
}

//---------------------------------------------------------------------------//
//                 end of tstParticle.cc
//---------------------------------------------------------------------------//
