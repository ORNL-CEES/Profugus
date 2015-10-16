//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstFission_Source.cc
 * \author Thomas M. Evans
 * \date   Tue May 06 11:54:26 2014
 * \brief  Fission_Source unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#include "comm/P_Stream.hh"
#include "geometry/Cartesian_Mesh.hh"
#include "../Fission_Source.hh"

#include "gtest/utils_gtest.hh"
#include "SourceTestBase.hh"

using def::X; using def::Y; using def::Z;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class FissionSourceTest : public SourceTestBase
{
    typedef SourceTestBase Base;

  protected:
    typedef profugus::Fission_Source         Fission_Source;
    typedef Fission_Source::SP_Fission_Sites SP_Fission_Sites;
    typedef Fission_Source::Fission_Site     Fission_Site;
    typedef std::shared_ptr<Fission_Source>  SP_Fission_Source;

    virtual int get_seed() const
    {
        return 3421;
    }

    // Initialize b_db
    virtual void init_db()
    {
        Base::init_db();

        if (nodes == 4)
            b_db->set("Np", 1250);
        else
            b_db->set("Np", 1000);
    }

    /*
      Test Lattice:

      |-------|-------|
      |       |       |
      |  UO2  |  H2O  |
      |       |       |
      |-------|-------|
      |       |       |
      |  H2O  |  UO2  |
      |       |       |
      |-------|-------|

      PIN: FUEL = 1
           MOD  = 0

      LATTICE:
           UO2 - 1
           H2O - 2
     */
    virtual void init_geometry()
    {
        typedef Geometry_t::SP_Array SP_Core;
        typedef Geometry_t::Array_t  Core_t;
        typedef Core_t::SP_Object    SP_Lattice;
        typedef Core_t::Object_t     Lattice_t;
        typedef Lattice_t::SP_Object SP_Pin_Cell;
        typedef Lattice_t::Object_t  Pin_Cell_t;

        // make pin cells
        SP_Pin_Cell uo2(std::make_shared<Pin_Cell_t>(1, 0.54, 0, 1.26, 14.28));
        SP_Pin_Cell h2o(std::make_shared<Pin_Cell_t>(0, 1.26, 14.28));

        // make lattice
        SP_Lattice lat(std::make_shared<Lattice_t>(2, 2, 1, 3));
        lat->assign_object(uo2, 1);
        lat->assign_object(h2o, 2);
        lat->id(0, 0, 0) = 2; // H20
        lat->id(1, 0, 0) = 1; // UO2
        lat->id(0, 1, 0) = 1; // UO2
        lat->id(1, 1, 0) = 2; // H2O

        // complete lattice
        lat->complete(0.0, 0.0, 0.0);

        // make the core
        SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));
        core->assign_object(lat, 0);
        core->complete(0.0, 0.0, 0.0);

        // make the b_geometry
        b_geometry = std::make_shared<Geometry_t>(core);
    }

    // 2 material definitions/1 group
    /*
     - Mat 0 -> Moderator
     - Mat 1 -> Fuel
     */
    virtual void init_physics()
    {
        RCP_XS xs(Teuchos::rcp(new XS_t()));
        xs->set(0, 1);

        XS_t::OneDArray total0(1, 1.1),   total1(1, 10.0);
        XS_t::TwoDArray scat0(1, 1, 0.9), scat1(1, 1, 2.1);

        XS_t::OneDArray chi1(1, 1.0);
        XS_t::OneDArray sigf1(1, 4.2);
        XS_t::OneDArray nusigf1(1, 2.4*4.2);

        xs->add(0, XS_t::TOTAL, total0);
        xs->add(0, 0, scat0);

        xs->add(1, XS_t::TOTAL, total1);
        xs->add(1, XS_t::CHI, chi1);
        xs->add(1, XS_t::NU_SIG_F, nusigf1);
        xs->add(1, XS_t::SIG_F, sigf1);

        XS_t::OneDArray bounds(b_group_bounds->group_bounds());

        xs->set_bounds(bounds);

        xs->complete();

        b_physics = std::make_shared<Physics_t>(b_db, xs);

        EXPECT_FALSE(b_physics->is_fissionable(0));
        EXPECT_TRUE(b_physics->is_fissionable(1));
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(FissionSourceTest, DefaultInitialize)
{
    EXPECT_FALSE(b_db->isParameter("init_fission_src"));

    // make fission source
    Fission_Source source(b_db, b_geometry, b_physics, b_rcon);

    EXPECT_DOUBLE_EQ(0., source.lower_coords()[X]);
    EXPECT_DOUBLE_EQ(0., source.lower_coords()[Y]);
    EXPECT_DOUBLE_EQ(0., source.lower_coords()[Z]);

    EXPECT_DOUBLE_EQ( 2.52, source.width()[X]);
    EXPECT_DOUBLE_EQ( 2.52, source.width()[Y]);
    EXPECT_DOUBLE_EQ(14.28, source.width()[Z]);

    // build the initial source
    EXPECT_TRUE(!profugus::Global_RNG::d_rng.assigned());
    source.build_initial_source();

    EXPECT_TRUE(!source.empty());

    if (nodes == 1)
    {
        EXPECT_EQ(1000, source.num_to_transport());
        EXPECT_EQ(1000, source.total_num_to_transport());
        EXPECT_EQ(1000, source.Np());
    }
    else if (nodes == 2)
    {
        EXPECT_EQ(500, source.num_to_transport());
        EXPECT_EQ(1000, source.total_num_to_transport());
        EXPECT_EQ(1000, source.Np());
    }
    else if (nodes == 4)
    {
        EXPECT_EQ(312, source.num_to_transport());
        EXPECT_EQ(1248, source.total_num_to_transport());
        EXPECT_EQ(1250, source.Np());
    }
}

//---------------------------------------------------------------------------//

TEST_F(FissionSourceTest, Initialize)
{
    Teuchos::Array<double> fbox(6, 0.0);
    fbox[1] = 2.52;
    fbox[3] = 2.52;
    fbox[5] = 14.28;
    b_db->set("init_fission_src", fbox);
    EXPECT_TRUE(b_db->isParameter("init_fission_src"));

    // make fission source
    Fission_Source source(b_db, b_geometry, b_physics, b_rcon);

    EXPECT_DOUBLE_EQ(0., source.lower_coords()[X]);
    EXPECT_DOUBLE_EQ(0., source.lower_coords()[Y]);
    EXPECT_DOUBLE_EQ(0., source.lower_coords()[Z]);

    EXPECT_DOUBLE_EQ( 2.52, source.width()[X]);
    EXPECT_DOUBLE_EQ( 2.52, source.width()[Y]);
    EXPECT_DOUBLE_EQ(14.28, source.width()[Z]);

    // build the initial source
    source.build_initial_source();

    EXPECT_TRUE(!source.empty());

    if (nodes == 1)
    {
        EXPECT_EQ(1000, source.num_to_transport());
        EXPECT_EQ(1000, source.total_num_to_transport());
        EXPECT_EQ(1000, source.Np());
    }
    else if (nodes == 2)
    {
        EXPECT_EQ(500, source.num_to_transport());
        EXPECT_EQ(1000, source.total_num_to_transport());
        EXPECT_EQ(1000, source.Np());
    }
    else if (nodes == 4)
    {
        EXPECT_EQ(312, source.num_to_transport());
        EXPECT_EQ(1248, source.total_num_to_transport());
        EXPECT_EQ(1250, source.Np());
    }
}

//---------------------------------------------------------------------------//

TEST_F(FissionSourceTest, Small)
{
    // set to a small number of particles and rebuild the source and test
    // particles
    b_db->set("Np", 12);

    // make fission source
    Fission_Source source(b_db, b_geometry, b_physics, b_rcon);
    source.build_initial_source();

    EXPECT_EQ(12, source.Np());
    EXPECT_EQ(12, source.total_num_to_transport());

    int ctr = 0;
    Fission_Source::SP_Particle p;
    while (!source.empty())
    {
        p = source.get_particle();
        ctr++;

        EXPECT_TRUE(p->alive());
        EXPECT_EQ(1, p->matid());
        EXPECT_EQ(1.0, p->wt());
    }

    if (nodes == 1)
    {
        EXPECT_EQ(12, ctr);
    }
    else if (nodes == 2)
    {
        EXPECT_EQ(6, ctr);
    }
    else if (nodes == 4)
    {
        EXPECT_EQ(3, ctr);
    }

    EXPECT_EQ(source.num_to_transport(), ctr);
    EXPECT_EQ(ctr, source.num_run());
    EXPECT_EQ(0, source.num_left());
}

//---------------------------------------------------------------------------//

TEST_F(FissionSourceTest, Initialize_Mesh)
{
    EXPECT_FALSE(b_db->isParameter("init_fission_src"));

    // make a mesh
    std::vector<double> r = {0.0, 1.26, 2.52};
    std::vector<double> z = {0.0, 5.0, 10.0, 14.28};
    auto mesh = std::make_shared<profugus::Cartesian_Mesh>(r, r, z);
    std::vector<double> rho = {0.0, 2.0, 1.5, 0.0,
                               0.0, 3.0, 2.5, 0.0,
                               0.0, 1.8, 1.2, 0.0};

    // make fission source
    Fission_Source source(b_db, b_geometry, b_physics, b_rcon);
    source.build_initial_source(mesh, rho);
    EXPECT_TRUE(!source.empty());

    if (nodes == 1)
    {
        EXPECT_EQ(999, source.num_to_transport());
        EXPECT_EQ(999, source.total_num_to_transport());
        EXPECT_EQ(1000, source.Np());
    }
    else if (nodes == 2)
    {
        EXPECT_TRUE(source.num_to_transport() >= 500 &&
                    source.num_to_transport() <= 501);
        EXPECT_EQ(1001, source.total_num_to_transport());
        EXPECT_EQ(1000, source.Np());
    }
    else if (nodes == 4)
    {
        EXPECT_TRUE(source.num_to_transport() >= 311 &&
                    source.num_to_transport() <= 314);
        EXPECT_EQ(1251, source.total_num_to_transport());
        EXPECT_EQ(1250, source.Np());
    }
}

//---------------------------------------------------------------------------//

TEST_F(FissionSourceTest, Small_Mesh)
{
    EXPECT_FALSE(b_db->isParameter("init_fission_src"));

    // make a mesh
    std::vector<double> r = {0.0, 1.26, 2.52};
    std::vector<double> z = {0.0, 5.0, 10.0, 14.28};
    auto mesh = std::make_shared<profugus::Cartesian_Mesh>(r, r, z);
    std::vector<double> rho = {0.0, 2.0, 1.5, 0.0,
                               0.0, 3.0, 2.5, 0.0,
                               0.0, 1.8, 1.2, 0.0};

    // set to a small number of particles and rebuild the source and test
    // particles
    b_db->set("Np", 12);

    // make fission source
    Fission_Source source(b_db, b_geometry, b_physics, b_rcon);
    source.build_initial_source(mesh, rho);
    EXPECT_TRUE(!source.empty());

    EXPECT_EQ(12, source.Np());

    int ctr = 0;
    Fission_Source::SP_Particle p;
    while (!source.empty())
    {
        p = source.get_particle();
        ctr++;

        EXPECT_TRUE(p->alive());
        EXPECT_EQ(1, p->matid());

        auto r = b_geometry->position(p->geo_state());

        if (r[1] < 1.26)
        {
            EXPECT_TRUE(r[0] > 1.26);
        }

        if (r[1] > 1.26)
        {
            EXPECT_TRUE(r[0] < 1.26);
        }

        if (nodes == 1)
        {
            EXPECT_EQ(1.0, p->wt());
        }
    }

    if (nodes == 1)
    {
        EXPECT_EQ(12, ctr);
    }

    EXPECT_EQ(source.num_to_transport(), ctr);
    EXPECT_EQ(ctr, source.num_run());
    EXPECT_EQ(0, source.num_left());
}

//---------------------------------------------------------------------------//

TEST_F(FissionSourceTest, FissionSrc)
{
    b_db->set("Np", 12);

    // make fission source
    Fission_Source source(b_db, b_geometry, b_physics, b_rcon);

    // build the initial source
    source.build_initial_source();
    EXPECT_EQ(nodes, source.num_streams());
    EXPECT_TRUE(profugus::Global_RNG::d_rng.assigned());
    EXPECT_TRUE(source.is_initial_source());

    // get a fission source container
    SP_Fission_Sites fsrc = source.create_fission_site_container();
    EXPECT_TRUE(static_cast<bool>(fsrc));
    EXPECT_TRUE(fsrc->empty());

    // add 3 sites to the container (6 total fissions)
    Fission_Site site;
    site.m = 1;
    site.r = Space_Vector(0.3, 2.1, 12.3);
    fsrc->push_back(site);
    fsrc->push_back(site);
    fsrc->push_back(site);
    site.m = 1;
    site.r = Space_Vector(0.25, 1.9, 12.1);
    fsrc->push_back(site);
    site.m = 1;
    site.r = Space_Vector(0.27, 1.8, 11.1);
    fsrc->push_back(site);
    fsrc->push_back(site);
    EXPECT_EQ(6, fsrc->size());
    double ref_hsrc = (-1.0 / 3.0 * (log10(1.0 / 3.0) / log10(2.0)))
                      + (-2.0 / 3.0 * (log10(2.0 / 3.0) / log10(2.0)));

    source.build_source(fsrc);
    EXPECT_FALSE(source.is_initial_source());
    EXPECT_TRUE(static_cast<bool>(fsrc));
    EXPECT_EQ(0, fsrc->size());
    EXPECT_TRUE(fsrc->empty());

    EXPECT_EQ(6, source.num_to_transport());
    EXPECT_EQ(6 * nodes, source.total_num_to_transport());
    EXPECT_EQ(2 * nodes, source.num_streams());

    int ctr       = 0;
    double tot_wt = 0.0;
    while (!source.empty())
    {
        Fission_Source::SP_Particle p = source.get_particle();
        ctr++;

        EXPECT_EQ(1, p->matid());
        EXPECT_TRUE(soft_equiv(p->wt(), 12.0 / (6 * nodes)));
        tot_wt += p->wt();

        EXPECT_EQ(p->rng().get_num(), profugus::Global_RNG::d_rng.get_num());

        if (nodes == 1)
        {
            EXPECT_EQ(1, p->rng().get_num());
        }
        else if (nodes == 2)
        {
            if (node == 0)
            {
                EXPECT_EQ(2, p->rng().get_num());
            }
            else
            {
                EXPECT_TRUE(p->rng().get_num() >= 3);
            }
        }
        else if (nodes == 4)
        {
            if (node == 0)
            {
                EXPECT_EQ(4, p->rng().get_num());
            }
            else if (node == 1)
            {
                EXPECT_TRUE(p->rng().get_num() >= 5);
            }
            else if (node == 2)
            {
                EXPECT_TRUE(p->rng().get_num() >= 6);
            }
            else if (node == 3)
            {
                EXPECT_TRUE(p->rng().get_num() >= 7);
            }
        }
    }

    profugus::global_sum(tot_wt);

    EXPECT_NEAR(12.0, tot_wt, 1.e-12);
    EXPECT_EQ(6, ctr);
    EXPECT_EQ(6, source.num_run());
    EXPECT_EQ(0, source.num_left());
}

//---------------------------------------------------------------------------//

TEST_F(FissionSourceTest, Distributions)
{
    if (nodes < 4)
         SKIP_TEST("Need at least 2 nodes; have " << nodes);

    b_db->set("Np", 50000);

    // make fission source
    Fission_Source source(b_db, b_geometry, b_physics, b_rcon);

    // build the initial source
    source.build_initial_source();
    EXPECT_EQ(nodes, source.num_streams());

    // bins
    std::vector<int> octant(8, 0);

    Fission_Source::SP_Particle p;
    while (!source.empty())
    {
        p = source.get_particle();

        EXPECT_TRUE(p->alive());
        EXPECT_EQ(1, p->matid());
        EXPECT_EQ(1.0, p->wt());

        Space_Vector omega = b_geometry->direction(p->geo_state());
        if (omega[0] > 0.0 && omega[1] > 0.0 && omega[2] > 0.0)
            octant[0]++;
        else if (omega[0] < 0.0 && omega[1] > 0.0 && omega[2] > 0.0)
            octant[1]++;
        else if (omega[0] > 0.0 && omega[1] < 0.0 && omega[2] > 0.0)
            octant[2]++;
        else if (omega[0] < 0.0 && omega[1] < 0.0 && omega[2] > 0.0)
            octant[3]++;
        else if (omega[0] > 0.0 && omega[1] > 0.0 && omega[2] < 0.0)
            octant[4]++;
        else if (omega[0] < 0.0 && omega[1] > 0.0 && omega[2] < 0.0)
            octant[5]++;
        else if (omega[0] > 0.0 && omega[1] < 0.0 && omega[2] < 0.0)
            octant[6]++;
        else if (omega[0] < 0.0 && omega[1] < 0.0 && omega[2] < 0.0)
            octant[7]++;
    }

    EXPECT_EQ(50000, source.total_num_to_transport());
    profugus::global_sum(&octant[0], 8);

    using profugus::setw; using profugus::fixed; using profugus::endl;
    using profugus::pcout;

        pcout << endl;
        pcout << "Angular Distributions of Fission Source (by octant)" << endl;
        for (int i = 0; i < 8; ++i)
        {
            double ans = static_cast<double>(octant[i]) / 50000.0;
            double ref = 1.0 / 8.0;
            double err = fabs(ans - ref) / ref;
            pcout << setw(3) << i << setw(12) << fixed << ans
                  << setw(12) << fixed << ref
                  << setw(12) << fixed << err
                  << endl;

            EXPECT_NEAR(ref, ans, .03);
        }
        pcout << endl;
}

//---------------------------------------------------------------------------//
//                 end of tstFission_Source.cc
//---------------------------------------------------------------------------//
