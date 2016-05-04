//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstParticle_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 15:22:15 2016
 * \brief  Test for Particle
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Particle_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"


//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(ParticleTest, random)
{
    Particle_Tester::test_randoms();
}

TEST(ParticleTest, groups)
{
    Particle_Tester::test_groups();
}

TEST(ParticleTest, matids)
{
    Particle_Tester::test_matids();
}

//---------------------------------------------------------------------------//
//                 end of tstParticle.cc
//---------------------------------------------------------------------------//
