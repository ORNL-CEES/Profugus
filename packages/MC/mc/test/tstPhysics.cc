//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstPhysics.cc
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
#include "../Sampler.hh"

using namespace std;

using def::X; using def::Y; using def::Z;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class PhysicsTest : public testing::Test
{
  protected:
    typedef profugus::RNG_Control             RNG_Control;
    typedef profugus::Physics                 Physics_t;
    typedef Physics_t::Geometry_t             Geometry_t;
    typedef Physics_t::SP_Geometry            SP_Geometry;
    typedef Physics_t::Particle_t             Particle;
    typedef Physics_t::Bank_t                 Bank_t;
    typedef Physics_t::SP_Particle            SP_Particle;
    typedef Physics_t::ParameterList_t        ParameterList_t;
    typedef Physics_t::RCP_Std_DB             RCP_Std_DB;
    typedef Physics_t::XS_t                   XS_t;
    typedef Physics_t::RCP_XS                 RCP_XS;
    typedef Physics_t::Fission_Site           Fission_Site;
    typedef Physics_t::Fission_Site_Container Fission_Site_Container;
    typedef shared_ptr<Physics_t>             SP_Physics;
    typedef Geometry_t::Space_Vector          Space_Vector;

  protected:

    void SetUp()
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

        build_xs();

        // make a rng
        int seed = 342412;
        profugus::RNG_Control control(seed);
        rng = control.rng();

        // make db
        db = Teuchos::rcp(new ParameterList_t("test"));
    }

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

        XS_t::OneDArray total(5);
        XS_t::TwoDArray scat(5, 5);

        double f[5] = {0.1, 0.4, 1.8, 5.7, 9.8};
        XS_t::OneDArray fission(begin(f), end(f));
        xs->add(1, XS_t::SIG_F, fission);

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

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(PhysicsTest, Collisions)
{
    // make a particle
    SP_Particle p(new Particle);
    geometry->initialize(
        Space_Vector(0.0, 0.0, 0.0), Space_Vector(1.0, 1.0, 1.0),
        p->geo_state());
    p->set_rng(rng);

    // check distributions from analog collisions
    {
        db->get("implicit_capture", false);

        // make a mg physics object
        SP_Physics physics(make_shared<Physics_t>(db, xs));
        physics->set_geometry(geometry);

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
            double rn = rng.ran();
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

                const Space_Vector &omega = geometry->direction(p->geo_state());

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
}

//---------------------------------------------------------------------------//
//                 end of tstPhysics.cc
//---------------------------------------------------------------------------//
