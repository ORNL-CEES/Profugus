//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_physics/test/tstMG_Physics.cc
 * \author Thomas M. Evans
 * \date   Tue Jan 04 13:22:55 2011
 * \brief  MG_Physics unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/nemesis_gtest.hh"

#include <vector>
#include <cmath>
#include <iomanip>
#include <cstring>

#include "utils/Definitions.hh"
#include "utils/SP.hh"
#include "database/Std_DB.hh"
#include "mc/RNG_Control.hh"
#include "mc/Sampler.hh"
#include "mc/Definitions.hh"
#include "geometry/rtk/RTK_Geometry.hh"

#include "../MG_Physics.hh"

using namespace std;

using mc::RNG_Control;
using mc::sampler::sample_discrete_CDF;

typedef shift::MG_Physics<denovo::RTK_Lattice>   MG_Physics_t;
typedef MG_Physics_t::Geometry_t                 Geometry_t;
typedef MG_Physics_t::SP_Geometry                SP_Geometry;
typedef MG_Physics_t::Particle_t                 Particle;
typedef MG_Physics_t::Bank_t                     Bank_t;
typedef MG_Physics_t::SP_Particle                SP_Particle;
typedef MG_Physics_t::Physics_State_t            State;
typedef MG_Physics_t::SP_Std_DB                  SP_Std_DB;
typedef MG_Physics_t::XS_DB_t                    XS_DB_t;
typedef MG_Physics_t::SP_XS_DB                   SP_XS_DB;
typedef MG_Physics_t::Fission_Site               Fission_Site;
typedef MG_Physics_t::Fission_Site_Container     Fission_Site_Container;
typedef nemesis::SP<MG_Physics_t>                SP_Physics;

typedef Geometry_t::Space_Vector Space_Vector;

using def::X; using def::Y; using def::Z;

int seed = 342412;

//---------------------------------------------------------------------------//
// HELPERS
//---------------------------------------------------------------------------//

SP_Geometry make_geometry()
{
    typedef Geometry_t::Array_t  Lattice_t;
    typedef Lattice_t::Object_t  Pin_Cell_t;
    typedef Geometry_t::SP_Array SP_Lattice;
    typedef Lattice_t::SP_Object SP_Pin_Cell;

    // make an infinite box
    SP_Pin_Cell box(new Pin_Cell_t(0, 100.0, 100.0));
    SP_Lattice  lat(new Lattice_t(1, 1, 1, 2));
    lat->id(0, 0, 0) = 1;
    lat->assign_object(box, 1);
    lat->complete(0.0, 0.0, 0.0);

    SP_Geometry geometry(new Geometry_t(lat));
    return geometry;
}

//---------------------------------------------------------------------------//
// TESTS
// ---------------------------------------------------------------------------//
// See scripts/tstMG_Physics.py for cross section data processing; the results
// are shown below:
/*
     | 1.2  0.0  0.0  0.0  0.0 |
     | 0.9  3.2  0.0  0.0  0.0 |
 S = | 0.4  2.8  6.9  1.5  0.0 |
     | 0.1  2.1  5.5  9.7  2.1 |
     | 0.0  0.2  1.3  6.6  9.9 |

   g    total      abs     scat        c
 ---------------------------------------
   0      5.2      2.6      2.6   0.5000
   1     11.4      3.1      8.3   0.7281
   2     18.2      4.5     13.7   0.7527
   3     29.9     12.1     17.8   0.5953
   4     27.3     15.3     12.0   0.4396

 Group-to-group scattering
 -------------------------
 0 -> 0 :   0.4615
 0 -> 1 :   0.3462
 0 -> 2 :   0.1538
 0 -> 3 :   0.0385
 0 -> 4 :   0.0000
 1 -> 0 :   0.0000
 1 -> 1 :   0.3855
 1 -> 2 :   0.3373
 1 -> 3 :   0.2530
 1 -> 4 :   0.0241
 2 -> 0 :   0.0000
 2 -> 1 :   0.0000
 2 -> 2 :   0.5036
 2 -> 3 :   0.4015
 2 -> 4 :   0.0949
 3 -> 0 :   0.0000
 3 -> 1 :   0.0000
 3 -> 2 :   0.0843
 3 -> 3 :   0.5449
 3 -> 4 :   0.3708
 4 -> 0 :   0.0000
 4 -> 1 :   0.0000
 4 -> 2 :   0.0000
 4 -> 3 :   0.1750
 4 -> 4 :   0.8250
 */

TEST(Physics, Collisions)
{
    // make a database
    SP_Std_DB db(new database::Std_DB("test"));
    db->new_key("num_groups", 5);
    db->new_key("Pn_order", 0);
    db->new_key("downscatter", false);

    // group boundaries
    vector<double> nbnd(6, 0.0);
    nbnd[0] = 100.0;
    nbnd[1] = 10.0;
    nbnd[2] = 1.0;
    nbnd[3] = 0.1;
    nbnd[4] = 0.01;
    nbnd[5] = 0.001;
    db->new_key("neutron_bnd", nbnd);

    // make cross sections
    SP_XS_DB xsdb(new XS_DB_t);
    {
        xsdb->set_num(1, 5);
        XS_DB_t::Vec_XS xs(5);
        for (int g = 0; g < 5; ++g)
            xs[g] = new XS_DB_t::XS(db, g);

        // totals
        xs[0]->sigma() = 5.2;
        xs[1]->sigma() = 11.4;
        xs[2]->sigma() = 18.2;
        xs[3]->sigma() = 29.9;
        xs[4]->sigma() = 27.3;

        // scattering
        xs[0]->sigma_s(0, 0) = 1.2;
        xs[1]->sigma_s(0, 0) = 0.9;
        xs[1]->sigma_s(1, 0) = 3.2;
        xs[2]->sigma_s(0, 0) = 0.4;
        xs[2]->sigma_s(1, 0) = 2.8;
        xs[2]->sigma_s(2, 0) = 6.9;
        xs[2]->sigma_s(3, 0) = 1.5;
        xs[3]->sigma_s(0, 0) = 0.1;
        xs[3]->sigma_s(1, 0) = 2.1;
        xs[3]->sigma_s(2, 0) = 5.5;
        xs[3]->sigma_s(3, 0) = 9.7;
        xs[3]->sigma_s(4, 0) = 2.1;
        xs[4]->sigma_s(1, 0) = 0.2;
        xs[4]->sigma_s(2, 0) = 1.3;
        xs[4]->sigma_s(3, 0) = 6.6;
        xs[4]->sigma_s(4, 0) = 9.9;

        // assign the cross sections
        for (int g = 0; g < 5; ++g)
        {
            xs[g]->complete();
            xsdb->assign(xs[g], 0);
        }
    }

    // make a geometry
    SP_Geometry geometry = make_geometry();

    // make a particle
    SP_Particle p(new Particle);
    geometry->change_direction(Space_Vector(1.0, 1.0, 1.0), p->geo_state());

    // make a rng
    RNG_Control control(seed);
    RNG_Control::RNG rng = control.rng();
    p->set_rng(rng);

    // check distributions from analog collisions
    {
        db->new_key("implicit_capture", false);

        // make a mg physics object
        SP_Physics physics(new MG_Physics_t(db, xsdb));
        physics->set_geometry(geometry);
        physics->set_rng(rng);

        // make a bank
        Bank_t bank(geometry, physics);

        int Np = 100000;

        int scat[5]   = {0};
        int abs[5]    = {0};
        int ng[5]     = {0};
        int octant[8] = {0};

        vector< vector<int> > g2g(5, vector<int>(5, 0));

        int ing, outg;

        // group CDF
        double cdf[] = {0.2, 0.4, 0.6, 0.8, 1.0};

        // physics state
        State &state = p->physics_state();

        for (int n = 0; n < Np; ++n)
        {
            // setup particle to do analog collision
            p->set_wt(1.0);
            p->set_event(mc::events::COLLISION);

            // sample a group
            double rn   = rng.ran();
            state.group = sample_discrete_CDF(5, cdf, rn);

            ing = state.group;
            ng[ing]++;

            // do the collision
            physics->collide(*p, bank);

            outg = state.group;

            if (p->event() == mc::events::ABSORPTION)
            {
                abs[ing]++;
            }
            else if (p->event() == mc::events::SCATTER)
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

    // cout << endl;
    // ut.passes("Analog collisions processed correctly.");

    // check distributions from implicit-capture collisions
    {
        db->set("implicit_capture", true);

        double c[] = {0.5000, 0.7281, 0.7527, 0.5953, 0.4396};

        // make a mg physics object
        SP_Physics physics(new MG_Physics_t(db, xsdb));

        physics->set_geometry(geometry);

        // physics state
        State &state = p->physics_state();

        // make a bank
        Bank_t bank(geometry, physics);

        int Np        = 100000;
        int octant[8] = {0};

        // group CDF
        double cdf[] = {0.2, 0.4, 0.6, 0.8, 1.0};

        for (int n = 0; n < Np; ++n)
        {

            // setup particle to do implicit capture
            p->set_wt(0.9);
            p->set_event(mc::events::COLLISION);

            // sample a group
            double rn   = rng.ran();
            state.group = sample_discrete_CDF(5, cdf, rn);
            int g       = state.group;

            // do the collision
            physics->collide(*p, bank);

            EXPECT_EQ(mc::events::IMPLICIT_CAPTURE, p->event());
            EXPECT_TRUE(soft_equiv(p->wt(), c[g] * 0.9, 1.0e-4));

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

    // cout << endl;
    // ut.passes("Implicit capture collisions processed correctly.");

    // cout << endl;
    // ut.passes("MG Physics collisions processed correctly.");
}

//---------------------------------------------------------------------------//

TEST(Data, Access)
{
    using mc::physics::TOTAL;
    using mc::physics::SCATTERING;
    using mc::physics::FISSION;

    // make a database
    SP_Std_DB db(new database::Std_DB("test"));
    db->new_key("num_groups", 5);
    db->new_key("Pn_order", 0);
    db->new_key("downscatter", false);

    // group boundaries
    vector<double> nbnd(6, 0.0);
    nbnd[0] = 100.0;
    nbnd[1] = 10.0;
    nbnd[2] = 1.0;
    nbnd[3] = 0.1;
    nbnd[4] = 0.01;
    nbnd[5] = 0.001;
    db->new_key("neutron_bnd", nbnd);

    // make cross sections
    SP_XS_DB xsdb(new XS_DB_t);
    {
        xsdb->set_num(3, 5);
        XS_DB_t::Vec_XS xs(5);
        for (int g = 0; g < 5; ++g)
            xs[g] = new XS_DB_t::XS(db, g);

        // totals
        xs[0]->sigma() = 5.2;
        xs[1]->sigma() = 11.4;
        xs[2]->sigma() = 18.2;
        xs[3]->sigma() = 29.9;
        xs[4]->sigma() = 27.3;

        // scattering
        xs[0]->sigma_s(0, 0) = 1.2;
        xs[1]->sigma_s(0, 0) = 0.9;
        xs[1]->sigma_s(1, 0) = 3.2;
        xs[2]->sigma_s(0, 0) = 0.4;
        xs[2]->sigma_s(1, 0) = 2.8;
        xs[2]->sigma_s(2, 0) = 6.9;
        xs[2]->sigma_s(3, 0) = 1.5;
        xs[3]->sigma_s(0, 0) = 0.1;
        xs[3]->sigma_s(1, 0) = 2.1;
        xs[3]->sigma_s(2, 0) = 5.5;
        xs[3]->sigma_s(3, 0) = 9.7;
        xs[3]->sigma_s(4, 0) = 2.1;
        xs[4]->sigma_s(1, 0) = 0.2;
        xs[4]->sigma_s(2, 0) = 1.3;
        xs[4]->sigma_s(3, 0) = 6.6;
        xs[4]->sigma_s(4, 0) = 9.9;

        // assign the cross sections
        for (int g = 0; g < 5; ++g)
        {
            xs[g]->complete();
            xsdb->assign(xs[g], 0);
        }
    }

    // make a state
    State state;
    int matid = 0;

    // test data access without fission
    {
        MG_Physics_t physics(db, xsdb);

        state.group = 0;
        EXPECT_SOFTEQ(5.2, physics.total(TOTAL, matid, state), 1.e-12);
        EXPECT_SOFTEQ(2.6, physics.total(SCATTERING, matid, state), 1.e-12);
        EXPECT_SOFTEQ(0.0, physics.total(FISSION, matid, state), 1.e-12);

        state.group = 1;
        EXPECT_SOFTEQ(11.4, physics.total(TOTAL, matid, state), 1.e-12);
        EXPECT_SOFTEQ(8.3, physics.total(SCATTERING, matid, state), 1.e-12);
        EXPECT_SOFTEQ(0.0, physics.total(FISSION, matid, state), 1.e-12);

        state.group = 2;
        EXPECT_SOFTEQ(18.2, physics.total(TOTAL, matid, state), 1.e-12);
        EXPECT_SOFTEQ(13.7, physics.total(SCATTERING, matid, state), 1.e-12);
        EXPECT_SOFTEQ(0.0, physics.total(FISSION, matid, state), 1.e-12);

        state.group = 3;
        EXPECT_SOFTEQ(29.9, physics.total(TOTAL, matid, state), 1.e-12);
        EXPECT_SOFTEQ(17.8, physics.total(SCATTERING, matid, state), 1.e-12);
        EXPECT_SOFTEQ(0.0, physics.total(FISSION, matid, state), 1.e-12);

        state.group = 4;
        EXPECT_SOFTEQ(27.3, physics.total(TOTAL, matid, state), 1.e-12);
        EXPECT_SOFTEQ(12.0, physics.total(SCATTERING, matid, state), 1.e-12);
        EXPECT_SOFTEQ(0.0, physics.total(FISSION, matid, state), 1.e-12);

    }

    // add a material with fission
    {
        XS_DB_t::Vec_XS xs(5);
        for (int g = 0; g < 5; ++g)
            xs[g] = new XS_DB_t::XS(db, g);

        double f[5] = {0.1, 0.4, 1.8, 5.7, 9.8};

        // totals
        xs[0]->sigma() = 5.2  + f[0];
        xs[1]->sigma() = 11.4 + f[1];
        xs[2]->sigma() = 18.2 + f[2];
        xs[3]->sigma() = 29.9 + f[3];
        xs[4]->sigma() = 27.3 + f[4];

        // scattering
        xs[0]->sigma_s(0, 0) = 1.2;
        xs[1]->sigma_s(0, 0) = 0.9;
        xs[1]->sigma_s(1, 0) = 3.2;
        xs[2]->sigma_s(0, 0) = 0.4;
        xs[2]->sigma_s(1, 0) = 2.8;
        xs[2]->sigma_s(2, 0) = 6.9;
        xs[2]->sigma_s(3, 0) = 1.5;
        xs[3]->sigma_s(0, 0) = 0.1;
        xs[3]->sigma_s(1, 0) = 2.1;
        xs[3]->sigma_s(2, 0) = 5.5;
        xs[3]->sigma_s(3, 0) = 9.7;
        xs[3]->sigma_s(4, 0) = 2.1;
        xs[4]->sigma_s(1, 0) = 0.2;
        xs[4]->sigma_s(2, 0) = 1.3;
        xs[4]->sigma_s(3, 0) = 6.6;
        xs[4]->sigma_s(4, 0) = 9.9;

        XS_DB_t::SP_FIS_XS fission(new XS_DB_t::FIS_XS(5));

        // assign the cross sections
        for (int g = 0; g < 5; ++g)
        {
            xs[g]->complete();
            xsdb->assign(xs[g], 1);

            fission->data(g, XS_DB_t::FIS_XS::SIGMA_F) = f[g];
        }
        xsdb->assign(fission, 1);
    }

    matid = 1;

    // test data access with fission
    {
        MG_Physics_t physics(db, xsdb);

        state.group = 0;
        EXPECT_SOFTEQ(5.3, physics.total(TOTAL, matid, state), 1.e-12);
        EXPECT_SOFTEQ(2.6, physics.total(SCATTERING, matid, state), 1.e-12);
        EXPECT_SOFTEQ(0.1, physics.total(FISSION, matid, state), 1.e-12);

        state.group = 1;
        EXPECT_SOFTEQ(11.8, physics.total(TOTAL, matid, state), 1.e-12);
        EXPECT_SOFTEQ(8.3, physics.total(SCATTERING, matid, state), 1.e-12);
        EXPECT_SOFTEQ(0.4, physics.total(FISSION, matid, state), 1.e-12);

        state.group = 2;
        EXPECT_SOFTEQ(20.0, physics.total(TOTAL, matid, state), 1.e-12);
        EXPECT_SOFTEQ(13.7, physics.total(SCATTERING, matid, state), 1.e-12);
        EXPECT_SOFTEQ(1.8, physics.total(FISSION, matid, state), 1.e-12);

        state.group = 3;
        EXPECT_SOFTEQ(35.6, physics.total(TOTAL, matid, state), 1.e-12);
        EXPECT_SOFTEQ(17.8, physics.total(SCATTERING, matid, state), 1.e-12);
        EXPECT_SOFTEQ(5.7, physics.total(FISSION, matid, state), 1.e-12);

        state.group = 4;
        EXPECT_SOFTEQ(37.1, physics.total(TOTAL, matid, state), 1.e-12);
        EXPECT_SOFTEQ(12.0, physics.total(SCATTERING, matid, state), 1.e-12);
        EXPECT_SOFTEQ(9.8, physics.total(FISSION, matid, state), 1.e-12);
    }


    // add a material with kappa data
    {
        XS_DB_t::Vec_XS xs(5);
        for (size_t g = 0; g < 5; ++g)
            xs[g] = new XS_DB_t::XS(db, g);

        XS_DB_t::SP_XS_Opt opt(new XS_DB_t::XS_Opt(5));

        for (size_t g = 0; g < 5; ++g)
        {
            xs[g]->sigma() = g + 0.5;

            for (size_t gp = 0; gp <= g; ++gp)
            {
                xs[g]->sigma_s(gp, 0) = 0.01 * g;
            }

            xs[g]->complete();
            xsdb->assign(xs[g], 2);

            opt->data(g, XS_DB_t::XS_Opt::KAPPA_SIGMA) = 100. * g + 1.;
        }

        xsdb->assign(opt, 2);
    }

    matid = 2;

    // test data access with fission
    {
        MG_Physics_t physics(db, xsdb);
        using mc::physics::KAPPA_SIGMA;

        state.group = 0;
        EXPECT_SOFTEQ(1.0, physics.total(KAPPA_SIGMA, matid, state), 1.e-12);

        state.group = 1;
        EXPECT_SOFTEQ(101.0, physics.total(KAPPA_SIGMA, matid, state), 1.e-12);

        state.group = 2;
        EXPECT_SOFTEQ(201.0, physics.total(KAPPA_SIGMA, matid, state), 1.e-12);

        state.group = 3;
        EXPECT_SOFTEQ(301.0, physics.total(KAPPA_SIGMA, matid, state), 1.e-12);

        state.group = 4;
        EXPECT_SOFTEQ(401.0, physics.total(KAPPA_SIGMA, matid, state), 1.e-12);
    }
}

//---------------------------------------------------------------------------//

TEST(Physics, initialization)
{
    // make a database
    SP_Std_DB db(new database::Std_DB("test"));
    db->new_key("Pn_order", 0);
    db->new_key("downscatter", false);

    // physics state
    State state;

    // 5 group neutron tests
    {
        db->new_key("num_groups", 5);

        // group boundaries
        vector<double> nbnd(6, 0.0);
        nbnd[0] = 100.0;
        nbnd[1] = 10.0;
        nbnd[2] = 1.0;
        nbnd[3] = 0.1;
        nbnd[4] = 0.01;
        nbnd[5] = 0.001;
        db->new_key("neutron_bnd", nbnd);

        // don't need to set cross sections to test initialization
        SP_XS_DB xsdb(new XS_DB_t);
        xsdb->set_num(1, 5);

        // make a mg physics object
        MG_Physics_t physics(db, xsdb);

        // check initializations
        physics.initialize(mc::NEUTRON, 100.0, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(0, state.group);

        physics.initialize(mc::NEUTRON, 99.99, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(0, state.group);

        physics.initialize(mc::NEUTRON, 10.01, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(0, state.group);

        physics.initialize(mc::NEUTRON, 10.0, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(0, state.group);

        physics.initialize(mc::NEUTRON, 9.99, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(1, state.group);

        physics.initialize(mc::NEUTRON, 1.01, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(1, state.group);

        physics.initialize(mc::NEUTRON, 1.0, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(1, state.group);

        physics.initialize(mc::NEUTRON, 0.99, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(2, state.group);

        physics.initialize(mc::NEUTRON, 0.101, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(2, state.group);

        physics.initialize(mc::NEUTRON, 0.1, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(2, state.group);

        physics.initialize(mc::NEUTRON, 0.099, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(3, state.group);

        physics.initialize(mc::NEUTRON, 0.011, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(3, state.group);

        physics.initialize(mc::NEUTRON, 0.01, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(3, state.group);

        physics.initialize(mc::NEUTRON, 0.0099, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(4, state.group);

        physics.initialize(mc::NEUTRON, 0.0011, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(4, state.group);

        physics.initialize(mc::NEUTRON, 0.001, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(4, state.group);

        // Check energy bounds
        EXPECT_EQ(0.001, physics.min_energy(mc::NEUTRON));
        EXPECT_EQ(100.0, physics.max_energy(mc::NEUTRON));

#ifdef CHECK_ON
        bool caught = false;
        try
        {
            physics.initialize(mc::NEUTRON, 0.00099, state);
        }
        catch (const nemesis::assertion &a)
        {
            caught = true;
        }
        EXPECT_TRUE(caught);
#endif
    }

    // 5 group photon tests
    {
        // group boundaries
        vector<double> pbnd(6, 0.0);
        pbnd[0] = 100.0;
        pbnd[1] = 10.0;
        pbnd[2] = 1.0;
        pbnd[3] = 0.1;
        pbnd[4] = 0.01;
        pbnd[5] = 0.001;
        db->new_key("photon_bnd", pbnd);
        db->erase< vector<double> >("neutron_bnd");

        // don't need to set cross sections to test initialization
        SP_XS_DB xsdb(new XS_DB_t);
        xsdb->set_num(1, 5);

        // make a mg physics object
        MG_Physics_t physics(db, xsdb);

        // check initializations

        physics.initialize(mc::PHOTON, 10.0, state);
        EXPECT_EQ(mc::PHOTON, physics.particle_type(state));
        EXPECT_EQ(0, state.group);

        physics.initialize(mc::PHOTON, 9.99, state);
        EXPECT_EQ(mc::PHOTON, physics.particle_type(state));
        EXPECT_EQ(1, state.group);

        physics.initialize(mc::PHOTON, 1.01, state);
        EXPECT_EQ(mc::PHOTON, physics.particle_type(state));
        EXPECT_EQ(1, state.group);

        physics.initialize(mc::PHOTON, 1.0, state);
        EXPECT_EQ(mc::PHOTON, physics.particle_type(state));
        EXPECT_EQ(1, state.group);

        // Check energy bounds
        EXPECT_EQ(0.001, physics.min_energy(mc::PHOTON));
        EXPECT_EQ(100.0, physics.max_energy(mc::PHOTON));

#ifdef CHECK_ON
        bool caught = false;
        try
        {
            physics.initialize(mc::PHOTON, 0.00099, state);
        }
        catch (const nemesis::assertion &a)
        {
            caught = true;
        }
        EXPECT_TRUE(caught);
#endif
    }

    // 5 group neutron/photon tests
    {
        // group boundaries
        vector<double> nbnd(4, 0.0);
        nbnd[0] = 20.0;
        nbnd[1] = 2.0;
        nbnd[2] = 0.02;
        nbnd[3] = 0.002;

        vector<double> pbnd(3, 0.0);
        pbnd[0] = 100.0;
        pbnd[1] = 1.0;
        pbnd[2] = 0.01;

        db->new_key("neutron_bnd", nbnd);
        db->set("photon_bnd", pbnd);

        // don't need to set cross sections to test initialization
        SP_XS_DB xsdb(new XS_DB_t);
        xsdb->set_num(1, 5);

        // make a mg physics object
        MG_Physics_t physics(db, xsdb);

        // check initializations

        physics.initialize(mc::NEUTRON, 20.0, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(0, state.group);

        physics.initialize(mc::NEUTRON, 10.0, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(0, state.group);

        physics.initialize(mc::NEUTRON, 2.0, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(0, state.group);

        physics.initialize(mc::NEUTRON, 0.021, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(1, state.group);

        physics.initialize(mc::NEUTRON, 0.002, state);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
        EXPECT_EQ(2, state.group);

        physics.initialize(mc::PHOTON, 9.1, state);
        EXPECT_EQ(mc::PHOTON, physics.particle_type(state));
        EXPECT_EQ(3, state.group);

        physics.initialize(mc::PHOTON, 1.01, state);
        EXPECT_EQ(mc::PHOTON, physics.particle_type(state));
        EXPECT_EQ(3, state.group);

        physics.initialize(mc::PHOTON, 1.0, state);
        EXPECT_EQ(mc::PHOTON, physics.particle_type(state));
        EXPECT_EQ(3, state.group);

        physics.initialize(mc::PHOTON,0.99, state);
        EXPECT_EQ(mc::PHOTON, physics.particle_type(state));
        EXPECT_EQ(4, state.group);

        physics.initialize(mc::PHOTON,0.01, state);
        EXPECT_EQ(mc::PHOTON, physics.particle_type(state));
        EXPECT_EQ(4, state.group);

        // Check energy bounds
        EXPECT_EQ(0.01, physics.min_energy(mc::PHOTON));
        EXPECT_EQ(100.0, physics.max_energy(mc::PHOTON));
        EXPECT_EQ(0.002, physics.min_energy(mc::NEUTRON));
        EXPECT_EQ(20.0, physics.max_energy(mc::NEUTRON));

#ifdef CHECK_ON
        bool caught = false;
        try
        {
            physics.initialize(mc::PHOTON, 0.00099, state);
        }
        catch (const nemesis::assertion &a)
        {
            caught = true;
        }
        EXPECT_TRUE(caught);
#endif
    }
}

//---------------------------------------------------------------------------//

TEST(State, Pack)
{
    // Testing MG Physics State pickle and restore

    // Make db and xsdb
    SP_Std_DB db(new database::Std_DB("test"));
    db->new_key("Pn_order", 0);
    db->new_key("downscatter", false);
    db->new_key("num_groups", 10);

    // Neutron group boundaries
    vector<double> nbnd(6, 0.0);
    nbnd[0] = 100.0;
    nbnd[1] = 10.0;
    nbnd[2] = 1.0;
    nbnd[3] = 0.1;
    nbnd[4] = 0.01;
    nbnd[5] = 0.001;
    db->new_key("neutron_bnd", nbnd);

    // Photon group boundaries
    vector<double> pbnd(6, 0.0);
    pbnd[0] = 200.0;
    pbnd[1] = 20.0;
    pbnd[2] = 2.0;
    pbnd[3] = 0.2;
    pbnd[4] = 0.02;
    pbnd[5] = 0.002;
    db->new_key("photon_bnd", pbnd);

    // Do not need to set cross sections to test packing/unpacking
    SP_XS_DB xsdb(new XS_DB_t);
    xsdb->set_num(1, 10);

    MG_Physics_t phys(db, xsdb);

    // Make a buffer
    vector<char> newbuffer;

    EXPECT_EQ(2 * sizeof(int), State::packed_bytes());

    State state;
    // Pack a state
    {
        phys.initialize(mc::PHOTON, 0.03, state);

        newbuffer.resize(state.packed_bytes());

        state.pack(&newbuffer[0]);
        EXPECT_EQ(8, state.group);
        EXPECT_EQ(0.02, phys.energy(state));
        EXPECT_EQ(mc::PHOTON, phys.particle_type(state));
    }

    // Change the state in the meantime
    State pstate;
    {
        phys.initialize(mc::NEUTRON, 1.01, pstate);
        EXPECT_EQ(1, pstate.group);
        EXPECT_EQ(1.0, phys.energy(pstate));
        EXPECT_EQ(mc::NEUTRON, phys.particle_type(pstate));

    }

    // Unpack the state
    State nstate;
    {
        nstate.unpack(&newbuffer[0]);

        EXPECT_EQ(8, nstate.group);
        EXPECT_EQ(0.02, phys.energy(nstate));
        EXPECT_EQ(mc::PHOTON, phys.particle_type(nstate));
    }

}

//---------------------------------------------------------------------------//

TEST(Fission, Sampling)
{
    // make a database
    SP_Std_DB db(new database::Std_DB("test"));
    db->new_key("num_groups", 5);
    db->new_key("Pn_order", 0);
    db->new_key("downscatter", false);

    // group boundaries
    vector<double> nbnd(6, 0.0);
    nbnd[0] = 100.0;
    nbnd[1] = 10.0;
    nbnd[2] = 1.0;
    nbnd[3] = 0.1;
    nbnd[4] = 0.01;
    nbnd[5] = 0.001;
    db->new_key("neutron_bnd", nbnd);

    // add a material with fission
    SP_XS_DB xsdb(new XS_DB_t);
    {
        xsdb->set_num(2, 5);
        XS_DB_t::Vec_XS xs(5);
        for (int g = 0; g < 5; ++g)
            xs[g] = new XS_DB_t::XS(db, g);

        double f[5] = {0.1, 0.4, 1.8, 5.7, 9.8};

        // totals
        xs[0]->sigma() = 5.2  + f[0];
        xs[1]->sigma() = 11.4 + f[1];
        xs[2]->sigma() = 18.2 + f[2];
        xs[3]->sigma() = 29.9 + f[3];
        xs[4]->sigma() = 27.3 + f[4];

        // scattering
        xs[0]->sigma_s(0, 0) = 1.2;
        xs[1]->sigma_s(0, 0) = 0.9;
        xs[1]->sigma_s(1, 0) = 3.2;
        xs[2]->sigma_s(0, 0) = 0.4;
        xs[2]->sigma_s(1, 0) = 2.8;
        xs[2]->sigma_s(2, 0) = 6.9;
        xs[2]->sigma_s(3, 0) = 1.5;
        xs[3]->sigma_s(0, 0) = 0.1;
        xs[3]->sigma_s(1, 0) = 2.1;
        xs[3]->sigma_s(2, 0) = 5.5;
        xs[3]->sigma_s(3, 0) = 9.7;
        xs[3]->sigma_s(4, 0) = 2.1;
        xs[4]->sigma_s(1, 0) = 0.2;
        xs[4]->sigma_s(2, 0) = 1.3;
        xs[4]->sigma_s(3, 0) = 6.6;
        xs[4]->sigma_s(4, 0) = 9.9;

        XS_DB_t::SP_FIS_XS fission(new XS_DB_t::FIS_XS(5));

        // assign the cross sections
        for (int g = 0; g < 5; ++g)
        {
            xs[g]->complete();
            xsdb->assign(xs[g], 1);

            fission->data(g, XS_DB_t::FIS_XS::SIGMA_F)    = f[g];
            fission->data(g, XS_DB_t::FIS_XS::NU_SIGMA_F) = 2.4 * f[g];
        }
        fission->data(0, XS_DB_t::FIS_XS::CHI) = 0.3770;
        fission->data(1, XS_DB_t::FIS_XS::CHI) = 0.4421;
        fission->data(2, XS_DB_t::FIS_XS::CHI) = 0.1809;

        xsdb->assign(fission, 1);
    }

    // make a geometry
    SP_Geometry geometry = make_geometry();

    MG_Physics_t physics(db, xsdb);
    physics.set_geometry(geometry);

    Fission_Site_Container fsites;

    // make a particle
    SP_Particle p(new Particle);
    geometry->initialize(Space_Vector(1.1, 0.5, 6.2),
                         Space_Vector(1.0, 1.0, 1.0),
                         p->geo_state());
    p->set_matid(1);
    p->set_wt(0.6);

    // make a rng
    RNG_Control control(seed);
    RNG_Control::RNG rng = control.rng();
    p->set_rng(rng);
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

    physics.set_rng(rng);

    // put particle in group 3
    p->physics_state().group = 3;
    p->physics_state().type  = mc::NEUTRON;

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
    p->physics_state().group = 4;
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
        State state;
        EXPECT_TRUE(physics.initialize_fission(1, rng, state));
        EXPECT_EQ(0, state.group);
        EXPECT_EQ(mc::NEUTRON, physics.particle_type(state));
    }

    // test initialization of particle from a fission site
    State s;
    EXPECT_TRUE(physics.initialize_fission(fsites[0], rng, s));
    EXPECT_EQ(1, s.group);
    EXPECT_EQ(mc::NEUTRON, s.type);

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
        Fission_Site_Container ufsc(4);
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
    EXPECT_FALSE(physics.initialize_fission(0, p->rng(), p->physics_state()));
    EXPECT_EQ(0, physics.sample_fission_site(*p, fsites, 0.1));
}

//---------------------------------------------------------------------------//
//                        end of tstMG_Physics.cc
//---------------------------------------------------------------------------//
