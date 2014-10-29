//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   acc/test/tstParticle.cc
 * \author Thomas M. Evans
 * \date   Tue Oct 28 13:43:27 2014
 * \brief  Test for Particle
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <memory>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "harness/DBC.hh"
#include "core/geometry/Mesh_Geometry.hh"
#include "core/mc/Uniform_Source.hh"
#include "core/mc/Box_Shape.hh"
#include "core/mc/Global_RNG.hh"
#include "core/mc/Group_Bounds.hh"
#include "core/mc/Physics.hh"

#include "../Particle.hh"
#include "ParticleTest.hh"


//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class ParticleTest : public testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::Physics                   Physics_t;
    typedef Physics_t::Geometry_t               Geometry_t;
    typedef profugus::Global_RNG::RNG_Control_t RNG_Control_t;

    typedef Physics_t::XS_t   XS_t;
    typedef Physics_t::RCP_XS RCP_XS;

    typedef Physics_t::ParameterList_t ParameterList_t;
    typedef Physics_t::RCP_Std_DB      RCP_Std_DB;

    typedef std::shared_ptr<Geometry_t>    SP_Geometry;
    typedef std::shared_ptr<Physics_t>     SP_Physics;
    typedef std::shared_ptr<RNG_Control_t> SP_RNG_Control;

    typedef std::shared_ptr<profugus::Group_Bounds> SP_Group_Bounds;

    typedef Physics_t::Space_Vector Space_Vector;

    typedef std::vector<double> Vec_Dbl;

    typedef profugus::Uniform_Source         Source;
    typedef Source::SP_Particle              SP_Particle;
    typedef std::shared_ptr<profugus::Shape> SP_Shape;

    typedef profugus::Mesh_Geometry       ACC_Geometry;
    typedef std::shared_ptr<ACC_Geometry> SP_ACC_Geometry;

  protected:
    void SetUp()
    {
        b_rcon = std::make_shared<RNG_Control_t>(this->get_seed());

        // set number of nodes
        node  = profugus::node();
        nodes = profugus::nodes();

        // set our other attributes
        init_group_bounds();
        init_db();
        init_geometry();
        init_physics();

        CHECK(b_physics);
        CHECK(b_geometry);
        b_physics->set_geometry(b_geometry);

        ENSURE(b_rcon);
        ENSURE(!b_db.is_null());
        ENSURE(b_geometry);
        ENSURE(b_physics);
        ENSURE(b_group_bounds);
    }

    void init_db()
    {
        REQUIRE(b_db.is_null());
        REQUIRE(b_group_bounds);

        b_db = Teuchos::rcp(new ParameterList_t("test"));
    }

    int get_seed() const
    {
        return 3421;
    }

    void init_group_bounds()
    {
        Vec_Dbl n_bounds;
        n_bounds.push_back(100.);
        n_bounds.push_back(0.001);

        b_group_bounds = std::make_shared<profugus::Group_Bounds>(n_bounds);
    }

    /*
      Test Geometry:

      |-------|-------|
      |       |       |
      |  H20  |  H2O  |
      |       |       |
      |-------|-------|
      |       |       |
      |  H2O  |  H20  |
      |       |       |
      |-------|-------|

     */
    void init_geometry()
    {
        typedef Geometry_t::SP_Array SP_Core;
        typedef Geometry_t::Array_t  Core_t;
        typedef Core_t::SP_Object    SP_Lattice;
        typedef Core_t::Object_t     Lattice_t;
        typedef Lattice_t::SP_Object SP_Pin_Cell;
        typedef Lattice_t::Object_t  Pin_Cell_t;

        // make cells
        SP_Pin_Cell h2o(std::make_shared<Pin_Cell_t>(3, 1.26, 14.28));

        // make lattice
        SP_Lattice lat(std::make_shared<Lattice_t>(2, 2, 1, 1));
        lat->assign_object(h2o, 0);

        // complete lattice
        lat->complete(0.0, 0.0, 0.0);

        // make the core
        SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));
        core->assign_object(lat, 0);
        core->complete(0.0, 0.0, 0.0);

        // make the b_geometry
        b_geometry = std::make_shared<Geometry_t>(core);

        // build a corresponding mesh geometry
        Vec_Dbl r = {0.0, 1.26, 2.52, 3.78};
        Vec_Dbl z = {0.0, 14.28};

        b_acc_geo = std::make_shared<ACC_Geometry>(r, r, z);
    }

    /*
     - Mat 3 -> H20
     */
    void init_physics()
    {
        RCP_XS xs(Teuchos::rcp(new XS_t()));
        xs->set(0, 1);

        XS_t::OneDArray total3(1, 1.1);
        XS_t::TwoDArray scat3(1, 1, 0.9);
        XS_t::OneDArray bounds(b_group_bounds->group_bounds());

        xs->add(3, XS_t::TOTAL, total3);
        xs->add(3, 0, scat3);
        xs->set_bounds(bounds);

        xs->complete();

        b_physics = std::make_shared<Physics_t>(b_db, xs);
    }

  protected:
    int             node;
    int             nodes;
    SP_RNG_Control  b_rcon;
    RCP_Std_DB      b_db;
    SP_Geometry     b_geometry;
    SP_ACC_Geometry b_acc_geo;
    SP_Physics      b_physics;
    SP_Group_Bounds b_group_bounds;

  protected:
    // >>> DATA
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(ParticleTest, regular_loop)
{
    b_db->set("Np", 48);

    // make a uniform source
    Source source(b_db, b_geometry, b_physics, b_rcon);

    // make a sampling shape (uniform)
    SP_Shape box(std::make_shared<profugus::Box_Shape>(
                     0.0, 2.52, 0.0, 2.52, 0.0, 14.28));

    // build the source
    source.build_source(box);

    EXPECT_TRUE(!source.empty());
    const int nodes = profugus::nodes();
    EXPECT_EQ(48 / nodes, source.num_to_transport());
    EXPECT_EQ(48, source.total_num_to_transport());
    EXPECT_EQ(48, source.Np());
    EXPECT_EQ(0, source.num_run());
    EXPECT_EQ(nodes, source.num_streams());

    int ctr = 0;
    while (!source.empty())
    {
        // get a particle
        SP_Particle p = source.get_particle();

        ctr++;

        EXPECT_TRUE(p->alive());
        EXPECT_EQ(1.0, p->wt());
        EXPECT_EQ(3, p->matid());
    }

    EXPECT_EQ(source.num_run(), ctr);
    EXPECT_EQ(source.num_to_transport(), ctr);
}

//---------------------------------------------------------------------------//

TEST_F(ParticleTest, load_particles)
{
    b_db->set("Np", 48);

    // make a uniform source
    Source source(b_db, b_geometry, b_physics, b_rcon);

    // make a sampling shape (uniform)
    SP_Shape box(std::make_shared<profugus::Box_Shape>(
                     0.0, 2.52, 0.0, 2.52, 0.0, 14.28));

    // build the source
    source.build_source(box);

    // get a vector of flattened particles
    acc::Vec_Particles particles;
    acc::set_size(particles, 10);

    // load the source
    acc::load_source(*acc_geometry, source, particles);

    // number of particles in grid
    EXPECT_EQ(60*128, particles.size());
    EXPECT_EQ(100, acc::rnd_numbers.size());

    for (int n = 0; n < 48; ++n)
    {
        // get a particle
        auto p = particles[n];

        EXPECT_TRUE(p.alive);
        EXPECT_EQ(1.0, p.wt);
        EXPECT_EQ(3, p.matid);
    }
}

//---------------------------------------------------------------------------//

TEST_F(ParticleTest, load_particles_gpu)
{
    b_db->set("Np", 48);

    // make a uniform source
    Source source(b_db, b_geometry, b_physics, b_rcon);

    // make a sampling shape (uniform)
    SP_Shape box(std::make_shared<profugus::Box_Shape>(
                     0.0, 2.52, 0.0, 2.52, 0.0, 14.28));

    // build the source
    source.build_source(box);

    // get a vector of flattened particles
    acc::Vec_Particles particles;
    acc::set_size(particles, 10);

    // load the source
    acc::load_source(*acc_geometry, source, particles);

    // number of particles in grid
    EXPECT_EQ(60*128, particles.size());
    EXPECT_EQ(100, acc::rnd_numbers.size());

    acc::Particle *p_ptr = &particles[0];

    loop_over_particles(p_ptr);
}

//---------------------------------------------------------------------------//
//                 end of tstParticle.cc
//---------------------------------------------------------------------------//
