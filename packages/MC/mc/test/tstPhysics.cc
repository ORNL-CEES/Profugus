//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/test/tstPhysics.cc
 * \author Thomas M. Evans
 * \date   Fri May 02 00:58:54 2014
 * \brief  Physics unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Physics.hh"

#include "gtest/utils_gtest.hh"

#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>

#include "Teuchos_RCP.hpp"
#include "utils/Definitions.hh"
#include "rng/RNG_Control.hh"
#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"
#include "../Sampler.hh"

using namespace std;

using def::X; using def::Y; using def::Z;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

template <class Geometry>
class PhysicsTest : public testing::Test
{
  protected:
    typedef Geometry                                    Geometry_t;
    typedef profugus::RNG_Control                       RNG_Control;
    typedef profugus::Physics<Geometry>                 Physics_t;
    typedef typename Physics_t::SP_Geometry             SP_Geometry;
    typedef typename Physics_t::Particle_t              Particle;
    typedef typename Physics_t::Bank_t                  Bank_t;
    typedef typename Physics_t::SP_Particle             SP_Particle;
    typedef typename Physics_t::ParameterList_t         ParameterList_t;
    typedef typename Physics_t::RCP_Std_DB              RCP_Std_DB;
    typedef typename Physics_t::XS_t                    XS_t;
    typedef typename Physics_t::RCP_XS                  RCP_XS;
    typedef typename Physics_t::Fission_Site            Fission_Site;
    typedef typename Physics_t::Fission_Site_Container  Fission_Site_Container;
    typedef shared_ptr<Physics_t>                       SP_Physics;
    typedef typename Geometry_t::Space_Vector           Space_Vector;

  protected:

    void SetUp()
    {
        build_geometry();

        build_xs();

        // make a rng
        int seed = 342412;
        profugus::RNG_Control control(seed);
        rng = control.rng();

        // make db
        db = Teuchos::rcp(new ParameterList_t("test"));
    }

    void build_geometry();

    void build_xs()
    {
        xs = Teuchos::rcp(new XS_t());
        xs->set(0, 5);

        vector<double> bnd(6, 0.0);
        bnd[0] = 100.0;
        bnd[1] = 10.0;
        bnd[2] = 1.0;
        bnd[3] = 0.1;
        bnd[4] = 0.01;
        bnd[5] = 0.001;
        xs->set_bounds(bnd);

        typename XS_t::OneDArray total(5);
        typename XS_t::TwoDArray scat(5, 5);

        double f[5] = {0.1, 0.4, 1.8, 5.7, 9.8};
        double c[5] = {0.3770, 0.4421, 0.1809, 0.0, 0.0};
        double n[5] = {2.4*f[0], 2.4*f[1], 2.4*f[2], 2.4*f[3], 2.4*f[4]};
        typename XS_t::OneDArray fission(begin(f), end(f));
        typename XS_t::OneDArray chi(begin(c), end(c));
        typename XS_t::OneDArray nus(begin(n), end(n));
        xs->add(1, XS_t::SIG_F, fission);
        xs->add(1, XS_t::NU_SIG_F, nus);
        xs->add(1, XS_t::CHI, chi);

        // mat 0
        total[0] = 5.2 ;
        total[1] = 11.4;
        total[2] = 18.2;
        total[3] = 29.9;
        total[4] = 27.3;
        xs->add(0, XS_t::TOTAL, total);

        // mat 1
        total[0] = 5.2  + f[0];
        total[1] = 11.4 + f[1];
        total[2] = 18.2 + f[2];
        total[3] = 29.9 + f[3];
        total[4] = 27.3 + f[4];
        xs->add(1, XS_t::TOTAL, total);

        scat(0, 0) = 1.2;
        scat(1, 0) = 0.9;
        scat(1, 1) = 3.2;
        scat(2, 0) = 0.4;
        scat(2, 1) = 2.8;
        scat(2, 2) = 6.9;
        scat(2, 3) = 1.5;
        scat(3, 0) = 0.1;
        scat(3, 1) = 2.1;
        scat(3, 2) = 5.5;
        scat(3, 3) = 9.7;
        scat(3, 4) = 2.1;
        scat(4, 1) = 0.2;
        scat(4, 2) = 1.3;
        scat(4, 3) = 6.6;
        scat(4, 4) = 9.9;
        xs->add(0, 0, scat);
        xs->add(1, 0, scat);

        xs->complete();
    }

  protected:
    SP_Geometry geometry;
    RCP_XS xs;
    RCP_Std_DB db;

    profugus::RNG_Control::RNG_t rng;
};

template <>
void PhysicsTest<profugus::Core>::build_geometry()
{
    typedef Geometry_t::Array_t  Core_t;
    typedef Core_t::Object_t     Lattice_t;
    typedef Lattice_t::Object_t  Pin_Cell_t;
    typedef Geometry_t::SP_Array SP_Core;
    typedef Core_t::SP_Object    SP_Lattice;
    typedef Lattice_t::SP_Object SP_Pin_Cell;

    // make an infinite box
    SP_Pin_Cell box(make_shared<Pin_Cell_t>(0, 100.0, 100.0));
    SP_Lattice  lat(make_shared<Lattice_t>(1, 1, 1, 2));
    lat->id(0, 0, 0) = 1;
    lat->assign_object(box, 1);
    lat->complete(0.0, 0.0, 0.0);

    SP_Core core(make_shared<Core_t>(1, 1, 1, 1));
    core->assign_object(lat, 0);
    core->complete(0.0, 0.0, 0.0);

    geometry = make_shared<Geometry_t>(core);
}

template <>
void PhysicsTest<profugus::Mesh_Geometry>::build_geometry()
{
    def::Vec_Dbl edges = {0.0, 100.0};
    auto matids = std::make_shared<def::Vec_Int>(def::Vec_Int({0}));

    geometry = std::make_shared<Geometry_t>(edges,edges,edges);
    geometry->set_matids(matids);
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

typedef ::testing::Types<profugus::Core,profugus::Mesh_Geometry> MyTypes;
TYPED_TEST_CASE(PhysicsTest, MyTypes);

TYPED_TEST(PhysicsTest, Collisions)
{
    typedef typename TestFixture::Particle      Particle;
    typedef typename TestFixture::SP_Particle   SP_Particle;
    typedef typename TestFixture::Physics_t     Physics_t;
    typedef typename TestFixture::SP_Physics    SP_Physics;
    typedef typename TestFixture::Space_Vector  Space_Vector;
    typedef typename TestFixture::Bank_t        Bank_t;

    // make a particle
    SP_Particle p(make_shared<Particle>());
    this->geometry->initialize(
        Space_Vector(50.0, 50.0, 50.0), Space_Vector(1.0, 1.0, 1.0),
        p->geo_state());
    p->set_rng(this->rng);
    p->set_matid(0);

    // check distributions from analog collisions
    {
        this->db->get("implicit_capture", false);

        // make a mg physics object
        SP_Physics physics(make_shared<Physics_t>(this->db, this->xs));
        physics->set_geometry(this->geometry);

        // make a bank
        Bank_t bank;

        int Np = 100000;

        int scat[5]   = {0};
        int abs[5]    = {0};
        int ng[5]     = {0};
        int octant[8] = {0};

        vector< vector<int> > g2g(5, vector<int>(5, 0));

        int ing, outg;

        // group CDF
        double cdf[] = {0.2, 0.4, 0.6, 0.8, 1.0};

        for (int n = 0; n < Np; ++n)
        {
            // setup particle to do analog collision
            p->set_wt(1.0);
            p->set_event(profugus::events::COLLISION);

            // sample a group
            double rn = this->rng.ran();
            int group = profugus::sampler::sample_discrete_CDF(5, cdf, rn);

            p->set_group(group);

            ing = group;
            ng[ing]++;

            // do the collision
            physics->collide(*p, bank);

            outg = p->group();

            if (p->event() == profugus::events::ABSORPTION)
            {
                abs[ing]++;
            }
            else if (p->event() == profugus::events::SCATTER)
            {
                scat[ing]++;
                g2g[ing][outg]++;

                const Space_Vector &omega =
                    this->geometry->direction(p->geo_state());

                // check octant distribution (isotropic scattering)
                if (omega[Z] > 0.0)
                {
                    if (omega[Y] > 0.0)
                    {
                        if (omega[X] > 0.0)
                        {
                            octant[0]++;
                        }
                        else
                        {
                            octant[1]++;
                        }
                    }
                    else
                    {
                        if (omega[X] > 0.0)
                        {
                            octant[2]++;
                        }
                        else
                        {
                            octant[3]++;
                        }
                    }
                }
                else
                {
                    if (omega[Y] > 0.0)
                    {
                        if (omega[X] > 0.0)
                        {
                            octant[4]++;
                        }
                        else
                        {
                            octant[5]++;
                        }
                    }
                    else
                    {
                        if (omega[X] > 0.0)
                        {
                            octant[6]++;
                        }
                        else
                        {
                            octant[7]++;
                        }
                    }
                }
            }
        }

        int tots         = 0, tota = 0;
        double ref[]     = {0.5000, 0.7281, 0.7527, 0.5953, 0.4396};
        double g2gr[][5] = {{0.4615, 0.3462, 0.1538, 0.0385, 0.0000},
                            {0.0000, 0.3855, 0.3373, 0.2530, 0.0241},
                            {0.0000, 0.0000, 0.5036, 0.4015, 0.0949},
                            {0.0000, 0.0000, 0.0843, 0.5449, 0.3708},
                            {0.0000, 0.0000, 0.0000, 0.1750, 0.8250}};

        cout.precision(4);
        cout << setw(3) << "g" << setw(10) << "abs" << setw(10) << "scat"
             << setw(10) << "c" << setw(10) << "diff" << endl;
        cout << "-------------------------------------------" << endl;
        for (int g = 0; g < 5; ++g)
        {
            double a = static_cast<double>(abs[g]) / ng[g];
            double s = static_cast<double>(scat[g]) / ng[g];

            cout << setw(3) << g
                 << setw(10) << fixed << a
                 << setw(10) << fixed << s
                 << setw(10) << fixed << ref[g]
                 << setw(10) << fixed << std::fabs(s - ref[g]) / ref[g]
                 << endl;

            EXPECT_SOFTEQ(s, ref[g], 0.008);

            tota += abs[g];
            tots += scat[g];
        }
        EXPECT_EQ(Np, tots + tota);

        cout << endl;
        cout << "Group-to-group Scattering" << endl;
        cout << "------------------------------------" << endl;
        for (int g = 0; g < 5; ++g)
        {
            for (int gp = 0; gp < 5; ++gp)
            {
                double s = static_cast<double>(g2g[g][gp]) / scat[g];
                double e = 0.0;
                if (g2gr[g][gp] > 0.0)
                    e = std::fabs(s - g2gr[g][gp]) / g2gr[g][gp];

                cout << setw(1) << g << " -> " << setw(1) << gp
                     << setw(10) << fixed << s
                     << setw(10) << fixed << g2gr[g][gp]
                     << setw(10) << fixed << e
                     << endl;

                EXPECT_SOFTEQ(s, g2gr[g][gp], 0.05);
            }
        }

        cout << endl;
        cout << "Isotropic scattering by octant" << endl;
        cout << "---------------------------------" << endl;
        int osum = 0;
        for (int i = 0; i < 8; ++i)
        {
            double o  = static_cast<double>(octant[i]) / tots;
            double r  = 1.0/8.0;
            osum     += octant[i];
            cout << setw(3) << i
                 << setw(10) << fixed << o
                 << setw(10) << fixed << r
                 << setw(10) << std::abs(o - r) / r << endl;

            EXPECT_SOFTEQ(o, r, 0.03);
        }
        EXPECT_EQ(tots, osum);
    }

    // check distributions from implicit-capture collisions
    {
        this->db->set("implicit_capture", true);

        double c[] = {0.5000, 0.7281, 0.7527, 0.5953, 0.4396};

        // make a mg physics object
        SP_Physics physics(make_shared<Physics_t>(this->db, this->xs));
        physics->set_geometry(this->geometry);

        // make a bank
        Bank_t bank;

        int Np        = 100000;
        int octant[8] = {0};

        // group CDF
        double cdf[] = {0.2, 0.4, 0.6, 0.8, 1.0};

        for (int n = 0; n < Np; ++n)
        {
            // setup particle to do implicit capture
            p->set_wt(0.9);
            p->set_event(profugus::events::COLLISION);

            // sample a group
            double rn   = this->rng.ran();
            int group   = profugus::sampler::sample_discrete_CDF(5, cdf, rn);

            p->set_group(group);

            int g = group;

            // do the collision
            physics->collide(*p, bank);

            EXPECT_EQ(profugus::events::IMPLICIT_CAPTURE, p->event());
            EXPECT_TRUE(soft_equiv(p->wt(), c[g] * 0.9, 1.0e-4));

            const Space_Vector &omega =
                this->geometry->direction(p->geo_state());

            // check octant distribution (isotropic scattering)
            if (omega[Z] > 0.0)
            {
                if (omega[Y] > 0.0)
                {
                    if (omega[X] > 0.0)
                    {
                        octant[0]++;
                    }
                    else
                    {
                        octant[1]++;
                    }
                }
                else
                {
                    if (omega[X] > 0.0)
                    {
                        octant[2]++;
                    }
                    else
                    {
                        octant[3]++;
                    }
                }
            }
            else
            {
                if (omega[Y] > 0.0)
                {
                    if (omega[X] > 0.0)
                    {
                        octant[4]++;
                    }
                    else
                    {
                        octant[5]++;
                    }
                }
                else
                {
                    if (omega[X] > 0.0)
                    {
                        octant[6]++;
                    }
                    else
                    {
                        octant[7]++;
                    }
                }
            }
        }

        cout << endl;
        cout << "Isotropic scattering by octant" << endl;
        cout << "---------------------------------" << endl;
        int osum = 0;
        for (int i = 0; i < 8; ++i)
        {
            double o  = static_cast<double>(octant[i]) / Np;
            double r  = 1.0/8.0;
            osum     += octant[i];
            cout << setw(3) << i
                 << setw(10) << fixed << o
                 << setw(10) << fixed << r
                 << setw(10) << std::abs(o - r) / r << endl;

            EXPECT_SOFTEQ(o, r, 0.03);
        }
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(PhysicsTest, Access)
{
    typedef typename TestFixture::Particle      Particle;
    typedef typename TestFixture::SP_Particle   SP_Particle;
    typedef typename TestFixture::Physics_t     Physics_t;
    typedef typename TestFixture::SP_Physics    SP_Physics;
    typedef typename TestFixture::Space_Vector  Space_Vector;

    using profugus::physics::TOTAL;
    using profugus::physics::SCATTERING;
    using profugus::physics::FISSION;

    // make a particle
    SP_Particle p(make_shared<Particle>());

    // physics
    Physics_t physics(this->db, this->xs);

    EXPECT_FALSE(physics.is_fissionable(0));
    EXPECT_TRUE(physics.is_fissionable(1));

    p->set_matid(0);

    // test data access without fission
    {

        p->set_group(0);
        EXPECT_SOFTEQ(5.2, physics.total(TOTAL, *p), 1.e-12);
        EXPECT_SOFTEQ(2.6, physics.total(SCATTERING, *p), 1.e-12);
        EXPECT_SOFTEQ(0.0, physics.total(FISSION, *p), 1.e-12);

        p->set_group(1);
        EXPECT_SOFTEQ(11.4, physics.total(TOTAL, *p), 1.e-12);
        EXPECT_SOFTEQ(8.3, physics.total(SCATTERING, *p), 1.e-12);
        EXPECT_SOFTEQ(0.0, physics.total(FISSION, *p), 1.e-12);

        p->set_group(2);
        EXPECT_SOFTEQ(18.2, physics.total(TOTAL, *p), 1.e-12);
        EXPECT_SOFTEQ(13.7, physics.total(SCATTERING, *p), 1.e-12);
        EXPECT_SOFTEQ(0.0, physics.total(FISSION, *p), 1.e-12);

        p->set_group(3);
        EXPECT_SOFTEQ(29.9, physics.total(TOTAL, *p), 1.e-12);
        EXPECT_SOFTEQ(17.8, physics.total(SCATTERING, *p), 1.e-12);
        EXPECT_SOFTEQ(0.0, physics.total(FISSION, *p), 1.e-12);

        p->set_group(4);
        EXPECT_SOFTEQ(27.3, physics.total(TOTAL, *p), 1.e-12);
        EXPECT_SOFTEQ(12.0, physics.total(SCATTERING, *p), 1.e-12);
        EXPECT_SOFTEQ(0.0, physics.total(FISSION, *p), 1.e-12);

    }

    p->set_matid(1);

    // test data access with fission
    {
        p->set_group(0);
        EXPECT_SOFTEQ(5.3, physics.total(TOTAL, *p), 1.e-12);
        EXPECT_SOFTEQ(2.6, physics.total(SCATTERING, *p), 1.e-12);
        EXPECT_SOFTEQ(0.1, physics.total(FISSION, *p), 1.e-12);

        p->set_group(1);
        EXPECT_SOFTEQ(11.8, physics.total(TOTAL, *p), 1.e-12);
        EXPECT_SOFTEQ(8.3, physics.total(SCATTERING, *p), 1.e-12);
        EXPECT_SOFTEQ(0.4, physics.total(FISSION, *p), 1.e-12);

        p->set_group(2);
        EXPECT_SOFTEQ(20.0, physics.total(TOTAL, *p), 1.e-12);
        EXPECT_SOFTEQ(13.7, physics.total(SCATTERING, *p), 1.e-12);
        EXPECT_SOFTEQ(1.8, physics.total(FISSION, *p), 1.e-12);

        p->set_group(3);
        EXPECT_SOFTEQ(35.6, physics.total(TOTAL, *p), 1.e-12);
        EXPECT_SOFTEQ(17.8, physics.total(SCATTERING, *p), 1.e-12);
        EXPECT_SOFTEQ(5.7, physics.total(FISSION, *p), 1.e-12);

        p->set_group(4);
        EXPECT_SOFTEQ(37.1, physics.total(TOTAL, *p), 1.e-12);
        EXPECT_SOFTEQ(12.0, physics.total(SCATTERING, *p), 1.e-12);
        EXPECT_SOFTEQ(9.8, physics.total(FISSION, *p), 1.e-12);
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(PhysicsTest, initialization)
{
    typedef typename TestFixture::Particle      Particle;
    typedef typename TestFixture::SP_Particle   SP_Particle;
    typedef typename TestFixture::Physics_t     Physics_t;
    typedef typename TestFixture::SP_Physics    SP_Physics;
    typedef typename TestFixture::Space_Vector  Space_Vector;

    // make a particle
    SP_Particle p(make_shared<Particle>());

    // physics
    Physics_t physics(this->db, this->xs);

    EXPECT_FALSE(physics.is_fissionable(0));
    EXPECT_TRUE(physics.is_fissionable(1));

    // check initializations
    physics.initialize(100.0, *p);
    EXPECT_EQ(0, p->group());

    physics.initialize(99.99, *p);
    EXPECT_EQ(0, p->group());

    physics.initialize(10.01, *p);
    EXPECT_EQ(0, p->group());

    physics.initialize(10.0, *p);
    EXPECT_EQ(0, p->group());

    physics.initialize(9.99, *p);
    EXPECT_EQ(1, p->group());

    physics.initialize(1.01, *p);
    EXPECT_EQ(1, p->group());

    physics.initialize(1.0, *p);
    EXPECT_EQ(1, p->group());

    physics.initialize(0.99, *p);
    EXPECT_EQ(2, p->group());

    physics.initialize(0.101, *p);
    EXPECT_EQ(2, p->group());

    physics.initialize(0.1, *p);
    EXPECT_EQ(2, p->group());

    physics.initialize(0.099, *p);
    EXPECT_EQ(3, p->group());

    physics.initialize(0.011, *p);
    EXPECT_EQ(3, p->group());

    physics.initialize(0.01, *p);
    EXPECT_EQ(3, p->group());

    physics.initialize(0.0099, *p);
    EXPECT_EQ(4, p->group());

    physics.initialize(0.0011, *p);
    EXPECT_EQ(4, p->group());

    physics.initialize(0.001, *p);
    EXPECT_EQ(4, p->group());

    // Check energy bounds
    EXPECT_EQ(0.001, physics.min_energy());
    EXPECT_EQ(100.0, physics.max_energy());

#ifdef CHECK_ON
    bool caught = false;
    try
    {
        physics.initialize(0.00099, *p);
    }
    catch (const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);
#endif
}



//---------------------------------------------------------------------------//

TYPED_TEST(PhysicsTest, fission_sampling)
{
    typedef typename TestFixture::Particle                  Particle;
    typedef typename TestFixture::SP_Particle               SP_Particle;
    typedef typename TestFixture::Physics_t                 Physics_t;
    typedef typename TestFixture::Space_Vector              Space_Vector;
    typedef typename TestFixture::Fission_Site_Container    FSC;

    Physics_t physics(this->db, this->xs);
    physics.set_geometry(this->geometry);

    FSC fsites;

    // make a particle
    SP_Particle p(make_shared<Particle>());
    this->geometry->initialize(Space_Vector(1.1, 0.5, 6.2),
                               Space_Vector(1.0, 1.0, 1.0),
                               p->geo_state());
    p->set_matid(1);
    p->set_wt(0.6);
    p->set_rng(this->rng);
    /*
     * First 10 random numbers in sequence:
     * 0.9709
     * 0.3771
     * 0.7536
     * 0.1897
     * 0.5297
     * 0.8803
     * 0.6286
     * 0.3288
     * 0.6362
     * 0.8904
     */

    // put particle in group 3
    p->set_group(3);

    // sampling fission with these setting should result in 1 fission event
    EXPECT_EQ(1, physics.sample_fission_site(*p, fsites, 1.04));
    EXPECT_EQ(1, fsites.size());
    EXPECT_EQ(1, fsites[0].m);
    EXPECT_EQ(1.1, fsites[0].r[0]);
    EXPECT_EQ(0.5, fsites[0].r[1]);
    EXPECT_EQ(6.2, fsites[0].r[2]);

    // this next one will fail
    p->geo_state().d_r = Space_Vector(1.2, 0.3, 6.6);
    EXPECT_EQ(0, physics.sample_fission_site(*p, fsites, 1.04));
    EXPECT_EQ(1, fsites.size());
    EXPECT_EQ(1, fsites[0].m);
    EXPECT_EQ(1.1, fsites[0].r[0]);
    EXPECT_EQ(0.5, fsites[0].r[1]);
    EXPECT_EQ(6.2, fsites[0].r[2]);

    // this one will pass
    p->set_group(4);
    p->set_wt(0.99);
    EXPECT_EQ(3, physics.sample_fission_site(*p, fsites, 0.2));
    EXPECT_EQ(4, fsites.size());
    EXPECT_EQ(1, fsites[0].m);
    EXPECT_EQ(1.1, fsites[0].r[0]);
    EXPECT_EQ(0.5, fsites[0].r[1]);
    EXPECT_EQ(6.2, fsites[0].r[2]);

    // there are 3 fission sites at this location
    EXPECT_EQ(1, fsites[1].m);
    EXPECT_EQ(1.2, fsites[1].r[0]);
    EXPECT_EQ(0.3, fsites[1].r[1]);
    EXPECT_EQ(6.6, fsites[1].r[2]);
    EXPECT_EQ(1, fsites[2].m);
    EXPECT_EQ(1.2, fsites[2].r[0]);
    EXPECT_EQ(0.3, fsites[2].r[1]);
    EXPECT_EQ(6.6, fsites[2].r[2]);
    EXPECT_EQ(1, fsites[2].m);
    EXPECT_EQ(1.2, fsites[2].r[0]);
    EXPECT_EQ(0.3, fsites[2].r[1]);
    EXPECT_EQ(6.6, fsites[2].r[2]);

    // test fission spectrum sample
    {
        EXPECT_TRUE(physics.initialize_fission(1, *p));
        EXPECT_EQ(0, p->group());
    }

    // test initialization of particle from a fission site
    EXPECT_TRUE(physics.initialize_fission(fsites[0], *p));
    EXPECT_EQ(1, p->group());

    EXPECT_EQ(4, fsites.size());

    EXPECT_EQ(1.1, physics.fission_site(fsites[0])[0]);
    EXPECT_EQ(0.5, physics.fission_site(fsites[0])[1]);
    EXPECT_EQ(6.2, physics.fission_site(fsites[0])[2]);
    EXPECT_EQ(1.2, physics.fission_site(fsites[1])[0]);
    EXPECT_EQ(0.3, physics.fission_site(fsites[1])[1]);
    EXPECT_EQ(6.6, physics.fission_site(fsites[1])[2]);

    cout << "\n Size of fission-site = " << sizeof(fsites[0]) << " bytes"
         << endl;

    EXPECT_EQ(sizeof(fsites[0]), physics.fission_site_bytes());

    // the size of the fission site container is 4 * 8 bytes (all of the
    // elements of the struct are aligned along 64-bit boundaries because of
    // the doubles in the space vector)
    EXPECT_EQ(32, physics.fission_site_bytes());

    // test packing of fission-container
    vector<char> packed;
    {
        packed.resize(physics.fission_site_bytes() * 4);
        memcpy(&packed[0], &fsites[0], packed.size());
    }

    {
        FSC ufsc(4);
        memcpy(&ufsc[0], &packed[0], packed.size());

        EXPECT_EQ(4, ufsc.size());
        EXPECT_EQ(1, ufsc[0].m);
        EXPECT_EQ(1.1, ufsc[0].r[0]);
        EXPECT_EQ(0.5, ufsc[0].r[1]);
        EXPECT_EQ(6.2, ufsc[0].r[2]);
        EXPECT_EQ(1, ufsc[1].m);
        EXPECT_EQ(1.2, ufsc[1].r[0]);
        EXPECT_EQ(0.3, ufsc[1].r[1]);
        EXPECT_EQ(6.6, ufsc[1].r[2]);
        EXPECT_EQ(1, ufsc[2].m);
        EXPECT_EQ(1.2, ufsc[2].r[0]);
        EXPECT_EQ(0.3, ufsc[2].r[1]);
        EXPECT_EQ(6.6, ufsc[2].r[2]);
        EXPECT_EQ(1, ufsc[3].m);
        EXPECT_EQ(1.2, ufsc[3].r[0]);
        EXPECT_EQ(0.3, ufsc[3].r[1]);
        EXPECT_EQ(6.6, ufsc[3].r[2]);
    }

    // test null ops for no fission
    p->set_matid(0);
    EXPECT_FALSE(physics.initialize_fission(0, *p));
    EXPECT_EQ(0, physics.sample_fission_site(*p, fsites, 0.1));
}

//---------------------------------------------------------------------------//
//                 end of tstPhysics.cc
//---------------------------------------------------------------------------//
