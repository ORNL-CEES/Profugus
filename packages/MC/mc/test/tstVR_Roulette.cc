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
#include "rng/RNG_Control.hh"
#include "../Definitions.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class VR_RouletteTest : public testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::Core                      Geometry_t;
    typedef profugus::VR_Roulette<Geometry_t>   Roulette;
    typedef Roulette::Physics_t                 Physics_t;
    typedef Roulette::Particle_t                Particle_t;
    typedef Roulette::Bank_t                    Bank_t;
    typedef Roulette::RCP_Std_DB                RCP_Std_DB;
    typedef Roulette::RNG_t                     RNG_t;
    typedef Roulette::SP_Geometry               SP_Geometry;
    typedef Roulette::SP_Physics                SP_Physics;

  protected:
    void SetUp()
    {
        db = Teuchos::rcp(new Roulette::ParameterList_t("test"));
        build_geometry();
        build_physics();

        seed = 23423;
        profugus::RNG_Control control(seed);

        auto ref = control.rng(12);
        rng      = control.rng(12);

        // reference random numbers
        refran.resize(10);
        for (auto &r : refran)
        {
            r = ref.ran();
        }

        p.set_rng(rng);
    }

    void build_geometry()
    {
        typedef Geometry_t::Array_t  Core_t;
        typedef Geometry_t::SP_Array SP_Core;
        typedef Core_t::Object_t     Lattice_t;
        typedef Core_t::SP_Object    SP_Lattice;
        typedef Lattice_t::Object_t  Pin_Cell_t;
        typedef Lattice_t::SP_Object SP_Pin_Cell;

        // make pin cells (mod_id, pitch, height)
        SP_Pin_Cell p1(std::make_shared<Pin_Cell_t>(0, 100., 100.));

        // make lattice (nx, ny, nz, num_objects)
        SP_Lattice lat(std::make_shared<Lattice_t>(1, 1, 1, 1));

        // assign pins
        lat->assign_object(p1, 0); // fuel pins

        // arrange pin-cells in lattice
        lat->id(0, 0, 0) = 0;

        // complete lattice
        lat->complete(0.0, 0.0, 0.0);

        // make core (nx, ny, nz, num_objects)
        SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));
        core->assign_object(lat, 0); // fuel pins
        core->id(0, 0, 0) = 0;
        core->complete(0.0, 0.0, 0.0);

        geometry = std::make_shared<Geometry_t>(core);
    }

    void build_physics()
    {
        typedef Physics_t::XS_t   XS_t;
        typedef Physics_t::RCP_XS RCP_XS;

        XS_t::OneDArray total(1, 0.9);
        XS_t::TwoDArray scat(1, 1, 0.8);

        double bnds[] = {200.0, 0.01};
        XS_t::OneDArray bounds(std::begin(bnds), std::end(bnds));

        RCP_XS xs(Teuchos::rcp(new XS_t));
        xs->set(0, 1);

        xs->add(3, XS_t::TOTAL, total);
        xs->add(3, 0, scat);
        xs->set_bounds(bounds);

        xs->complete();

        physics = std::make_shared<Physics_t>(db, xs);
        physics->set_geometry(geometry);
    }

  protected:
    // >>> DATA

    RCP_Std_DB  db;
    SP_Geometry geometry;
    SP_Physics  physics;

    int seed;

    Particle_t p;
    Bank_t     b;
    RNG_t      rng;

    std::vector<double> refran;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(VR_RouletteTest, default_settings)
{
    Roulette vr(db);
    EXPECT_EQ(0.5,  vr.weight_survival());
    EXPECT_EQ(0.25, vr.weight_cutoff());

    p.live();
    p.set_event(profugus::events::COLLISION);

    // above cutoff
    p.set_wt(0.2501);
    vr.post_collision(p, b);
    EXPECT_TRUE(p.alive());
    EXPECT_EQ(0.2501, p.wt());
    EXPECT_EQ(profugus::events::COLLISION, p.event());

    // below cutoff (particle rouletted)
    p.set_wt(0.2499);
    vr.post_collision(p, b);
    EXPECT_FALSE(p.alive());
    EXPECT_EQ(0.2499, p.wt());
    EXPECT_EQ(profugus::events::ROULETTE_KILLED, p.event());

    p.live();

    // below cutoff (particle survives)
    p.set_wt(0.058);
    vr.post_collision(p, b);
    EXPECT_TRUE(p.alive());
    EXPECT_EQ(0.5, p.wt());
    EXPECT_EQ(profugus::events::ROULETTE_SURVIVE, p.event());

    // below cutoff (particle rouletted)
    p.set_wt(0.20);
    vr.post_collision(p, b);
    EXPECT_FALSE(p.alive());
    EXPECT_EQ(0.20, p.wt());
    EXPECT_EQ(profugus::events::ROULETTE_KILLED, p.event());

    p.live();

    // below cutoff (particle rouletted)
    p.set_wt(0.212501);
    vr.post_collision(p, b);
    EXPECT_FALSE(p.alive());
    EXPECT_EQ(0.212501, p.wt());
    EXPECT_EQ(profugus::events::ROULETTE_KILLED, p.event());

    EXPECT_EQ(0, b.num_particles());
}

//---------------------------------------------------------------------------//

TEST_F(VR_RouletteTest, zero_cutoff)
{
    db->set("weight_cutoff", 0.0);
    Roulette vr(db);
    EXPECT_EQ(0.0, vr.weight_survival());
    EXPECT_EQ(0.0, vr.weight_cutoff());

    p.live();
    p.set_event(profugus::events::COLLISION);

    // above cutoff
    p.set_wt(0.2501);
    vr.post_collision(p, b);
    EXPECT_TRUE(p.alive());
    EXPECT_EQ(0.2501, p.wt());
    EXPECT_EQ(profugus::events::COLLISION, p.event());

    // above cutoff
    p.set_wt(0.000001);
    vr.post_collision(p, b);
    EXPECT_TRUE(p.alive());
    EXPECT_EQ(0.000001, p.wt());
    EXPECT_EQ(profugus::events::COLLISION, p.event());

    // should be no particles in the bank
    EXPECT_EQ(0, b.num_particles());
}

//---------------------------------------------------------------------------//

TEST_F(VR_RouletteTest, assert)
{
    db->set("weight_cutoff", 0.1);
    db->set("weight_survival", 0.09);

    EXPECT_THROW(Roulette vr(db), profugus::assertion);
}

//---------------------------------------------------------------------------//

TEST_F(VR_RouletteTest, construct1)
{
    db->set("weight_cutoff", 0.0);
    db->set("weight_survival", 0.09);
    Roulette vr(db);
    EXPECT_EQ(0.09, vr.weight_survival());
    EXPECT_EQ(0.0,  vr.weight_cutoff());
}

//---------------------------------------------------------------------------//

TEST_F(VR_RouletteTest, construct2)
{
    db->set("weight_cutoff", 0.01);
    Roulette vr(db);
    EXPECT_EQ(0.02,  vr.weight_survival());
    EXPECT_EQ(0.01,  vr.weight_cutoff());
}

//---------------------------------------------------------------------------//

TEST_F(VR_RouletteTest, construct3)
{
    db->set("weight_survival", 0.9);
    Roulette vr(db);
    EXPECT_EQ(0.9,  vr.weight_survival());
    EXPECT_EQ(0.25, vr.weight_cutoff());
}

//---------------------------------------------------------------------------//
//                 end of tstVR_Roulette.cc
//---------------------------------------------------------------------------//
