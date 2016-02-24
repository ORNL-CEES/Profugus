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
// Test fixture
//---------------------------------------------------------------------------//
class ParticleTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS

  protected:
    void SetUp()
    {
        /* * */
    }

  protected:
    // >>> DATA
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(ParticleTest, random)
{
    int num_vals = 16;
    std::vector<double> rands(num_vals,0.0);

    cuda_mc::Particle_Tester::test_randoms(rands);

    for( auto val : rands )
    {
        EXPECT_TRUE( val > 0.0 );
        EXPECT_TRUE( val < 1.0 );
    }
}

TEST_F(ParticleTest, groups)
{
    int num_vals = 16;
    std::vector<int> groups_in(num_vals,0.0);
    for( int i = 0; i < num_vals; ++i )
        groups_in[i] = i % 4;
    std::vector<int> groups_out(num_vals,0.0);

    cuda_mc::Particle_Tester::test_groups(groups_in,groups_out);

    for( int i = 0; i < num_vals; ++i )
    {
        EXPECT_EQ( i % 4, groups_out[i] );
    }
}

TEST_F(ParticleTest, matids)
{
    int num_vals = 16;
    std::vector<int> matids_in(num_vals,0.0);
    for( int i = 0; i < num_vals; ++i )
        matids_in[i] = i % 4;
    std::vector<int> matids_out(num_vals,0.0);

    cuda_mc::Particle_Tester::test_matids(matids_in,matids_out);

    for( int i = 0; i < num_vals; ++i )
    {
        EXPECT_EQ( i % 4, matids_out[i] );
    }
}

//---------------------------------------------------------------------------//
//                 end of tstParticle.cc
//---------------------------------------------------------------------------//
