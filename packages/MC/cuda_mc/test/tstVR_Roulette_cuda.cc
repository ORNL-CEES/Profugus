//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstVR_Roulette.cc
 * \author Thomas M. Evans
 * \date   Fri May 09 13:29:04 2014
 * \brief  VR_Roulette unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../VR_Roulette.hh"

#include "gtest/utils_gtest.hh"

#include <vector>
#include <memory>
#include "../Definitions.hh"
#include "rng/RNG_Control.hh"
#include "Particle_Vector_Tester.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class VR_RouletteTest : public testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef cuda_profugus::Mesh_Geometry           Geometry_t;
    typedef cuda_profugus::VR_Roulette<Geometry_t> Roulette;
    typedef Roulette::Bank_t                       Bank_t;
    typedef profugus::RNG_Control                  RNG_Control;

  protected:
    void SetUp()
    {
        int seed = 342412;
        profugus::RNG_Control control(seed);
        rng = control.rng();
    }

  protected:
    Teuchos::ParameterList db;
    profugus::RNG_Control::RNG_t rng;
    cuda::Shared_Device_Ptr<Bank_t> bank;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(VR_RouletteTest, default_settings)
{
    int num_particle = 1;
    Particle_Vector_Tester tester( num_particle, rng );

    tester.live();
    Teuchos::Array<cuda_profugus::events::Event> events( 
	num_particle, cuda_profugus::events::COLLISION );
    tester.set_event(events);
    tester.sort_by_event();

    Roulette vr(db);
    EXPECT_EQ(0.5,  vr.weight_survival());
    EXPECT_EQ(0.25, vr.weight_cutoff());

    // above cutoff. all should live
    tester.set_wt(0.2501);
    vr.post_collision(tester.get_vector(), bank);
    auto alive = tester.alive();
    for ( auto i : alive ) EXPECT_TRUE( i );
    auto wt = tester.wt();
    for ( auto i : wt ) EXPECT_EQ( 0.2501, i );
    auto event = tester.event();
    for ( auto i : event ) EXPECT_EQ( cuda_profugus::events::COLLISION, i );

    // below cutoff. will survive
    tester.set_wt(0.2499);
    vr.post_collision(tester.get_vector(), bank);
    alive = tester.alive();
    for ( auto i : alive ) EXPECT_TRUE( i );
    wt = tester.wt();
    for ( auto i : wt ) EXPECT_EQ( 0.5, i );
    event = tester.event();
    for ( auto i : event ) EXPECT_EQ( cuda_profugus::events::ROULETTE_SURVIVE, i );

    tester.live();
    tester.set_event(events);
    tester.sort_by_event();

    // below cutoff. will die
    tester.set_wt(0.018);
    vr.post_collision(tester.get_vector(), bank);
    alive = tester.alive();
    for ( auto i : alive ) EXPECT_FALSE( i );
    wt = tester.wt();
    for ( auto i : wt ) EXPECT_EQ( 0.018, i );
    event = tester.event();
    for ( auto i : event ) EXPECT_EQ( cuda_profugus::events::ROULETTE_KILLED, i );

    tester.live();
    tester.set_event(events);
    tester.sort_by_event();

    // below cutoff. some should roulette
    tester.set_wt(0.20);
    vr.post_collision(tester.get_vector(), bank);
    alive = tester.alive();
    for ( auto i : alive ) EXPECT_FALSE( i );
    wt = tester.wt();
    for ( auto i : wt ) EXPECT_EQ( 0.2, i );
    event = tester.event();
    for ( auto i : event ) EXPECT_EQ( cuda_profugus::events::ROULETTE_KILLED, i );

    tester.live();
    tester.set_event(events);
    tester.sort_by_event();

    // below cutoff. will survive
    tester.set_wt(0.212501);
    vr.post_collision(tester.get_vector(), bank);
    alive = tester.alive();
    for ( auto i : alive ) EXPECT_TRUE( i );
    wt = tester.wt();
    for ( auto i : wt ) EXPECT_EQ( 0.5, i );
    event = tester.event();
    for ( auto i : event ) EXPECT_EQ( cuda_profugus::events::ROULETTE_SURVIVE, i );
}

//---------------------------------------------------------------------------//

TEST_F(VR_RouletteTest, zero_cutoff)
{
    int num_particle = 10;
    Particle_Vector_Tester tester( num_particle, rng );

    tester.live();
    Teuchos::Array<cuda_profugus::events::Event> events( 
	num_particle, cuda_profugus::events::COLLISION );
    tester.set_event(events);
    tester.sort_by_event();

    db.set("weight_cutoff", 0.0);
    Roulette vr(db);
    EXPECT_EQ(0.0, vr.weight_survival());
    EXPECT_EQ(0.0, vr.weight_cutoff());

    // above cutoff
    tester.set_wt(0.2501);
    vr.post_collision(tester.get_vector(), bank);
    auto alive = tester.alive();
    for ( auto i : alive ) EXPECT_TRUE( i );
    auto wt = tester.wt();
    for ( auto i : wt ) EXPECT_EQ( 0.2501, i );
    auto event = tester.event();
    for ( auto i : event ) EXPECT_EQ( cuda_profugus::events::COLLISION, i );

    tester.live();
    tester.set_event(events);
    tester.sort_by_event();

    // above cutoff
    tester.set_wt(0.000001);
    vr.post_collision(tester.get_vector(), bank);
    alive = tester.alive();
    for ( auto i : alive ) EXPECT_TRUE( i );
    wt = tester.wt();
    for ( auto i : wt ) EXPECT_EQ( 0.000001, i );
    event = tester.event();
    for ( auto i : event ) EXPECT_EQ( cuda_profugus::events::COLLISION, i );

    tester.live();
    tester.set_event(events);
    tester.sort_by_event();
}

//---------------------------------------------------------------------------//

TEST_F(VR_RouletteTest, construct1)
{
    db.set("weight_cutoff", 0.0);
    db.set("weight_survival", 0.09);
    Roulette vr(db);
    EXPECT_EQ(0.09, vr.weight_survival());
    EXPECT_EQ(0.0,  vr.weight_cutoff());
}

//---------------------------------------------------------------------------//

TEST_F(VR_RouletteTest, construct2)
{
    db.set("weight_cutoff", 0.01);
    Roulette vr(db);
    EXPECT_EQ(0.02,  vr.weight_survival());
    EXPECT_EQ(0.01,  vr.weight_cutoff());
}

//---------------------------------------------------------------------------//

TEST_F(VR_RouletteTest, construct3)
{
    db.set("weight_survival", 0.9);
    Roulette vr(db);
    EXPECT_EQ(0.9,  vr.weight_survival());
    EXPECT_EQ(0.25, vr.weight_cutoff());
}

//---------------------------------------------------------------------------//
//                 end of tstVR_Roulette.cc
//---------------------------------------------------------------------------//
