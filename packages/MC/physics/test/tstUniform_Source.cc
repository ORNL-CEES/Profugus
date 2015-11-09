//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstUniform_Source.cc
 * \author Thomas M. Evans
 * \date   Wed May 07 15:24:31 2014
 * \brief  Uniform_Source test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>

#include "../Uniform_Source.hh"
#include "../Box_Shape.hh"

#include "gtest/utils_gtest.hh"

#include "SourceTestBase.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class UniformSourceTest : public SourceTestBase
{
    typedef SourceTestBase Base;

  protected:
    typedef profugus::Uniform_Source         Source;
    typedef Source::Particle_t              Particle_t;
    typedef std::shared_ptr<profugus::Shape> SP_Shape;

    virtual int get_seed() const
    {
        return 3421;
    }

    virtual void init_group_bounds()
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
    virtual void init_geometry()
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
    }

    /*
     - Mat 3 -> H20
     */
    virtual void init_physics()
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
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(UniformSourceTest, build_and_run)
{
    int np = 48;
    b_db->set("Np", np);

    // make a sampling shape (uniform)
    SP_Shape box(std::make_shared<profugus::Box_Shape>(
                     0.0, 2.52, 0.0, 2.52, 0.0, 14.28));

    // make a uniform source
    Source source(b_db, b_geometry, b_physics, box);

    const int nodes = profugus::nodes();
    EXPECT_EQ(48 / nodes, source.num_to_transport());
    EXPECT_EQ(48, source.total_num_to_transport());
    EXPECT_EQ(48, source.Np());

    int ctr = 0;
    for ( int i = 0; i < np; ++i )
    {
        // get a particle
        Particle_t p = source.get_particle(i);

        ctr++;

        EXPECT_TRUE(p.alive());
        EXPECT_EQ(1.0, p.wt());
        EXPECT_EQ(3, p.matid());
    }

    EXPECT_EQ(source.num_to_transport(), ctr);
}

//---------------------------------------------------------------------------//
//                 end of tstUniform_Source.cc
//---------------------------------------------------------------------------//
