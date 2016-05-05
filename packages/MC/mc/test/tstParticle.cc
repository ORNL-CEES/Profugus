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
#include "geometry/RTK_Geometry.hh"

typedef profugus::Particle<profugus::Core> Particle;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Particle, construction)
{
    Particle p;

    // access the particle metadata
    EXPECT_EQ(0, p.metadata().size());
}

//---------------------------------------------------------------------------//

TEST(Particle, metadata)
{
    // add some metadata
    auto fm_cell = Particle::Metadata::new_pod_member<int>("fm_birth_cell");

    // make the particle
    {
        Particle p, q;

        EXPECT_EQ(1, p.metadata().size());

        auto &metadata = p.metadata();
        metadata.access<int>(fm_cell) = 101;

        EXPECT_EQ(101, p.metadata().access<int>(fm_cell));
    }

    Particle::Metadata::reset();

    Particle p;
    EXPECT_EQ(0, p.metadata().size());
}

//---------------------------------------------------------------------------//
//                 end of tstParticle.cc
//---------------------------------------------------------------------------//
