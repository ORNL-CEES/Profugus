//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstDomain_Transporter.cc
 * \author Thomas M. Evans
 * \date   Mon May 12 14:14:23 2014
 * \brief  Domain_Transporter unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <sstream>
#include <cmath>

#include "../Domain_Transporter.hh"
#include "../VR_Roulette.hh"

#include "gtest/utils_gtest.hh"

#include "TransporterTestBase.hh"

using namespace std;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class Domain_TransporterTest : public TransporterTestBase
{
  protected:
    typedef profugus::Domain_Transporter Transporter;

  protected:

    void init_vr()
    {
        db->set("weight_cutoff", 0.001);
        var_red = std::make_shared<profugus::VR_Roulette>(db);
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Domain_TransporterTest, transport)
{
    Transporter transporter;
    transporter.set(geometry, physics);
    transporter.set(var_red);
    transporter.set(tallier);

    // make a particle and bank
    Particle_t p;
    Bank_t     bank;

    // make random number
    auto rng = rcon->rng(15);
    p.set_rng(rng);

    // geometry and physics state
    Particle_t::Geo_State_t &geo = p.geo_state();

    // events
    int esc = 0, rk = 0;

    // number of particles
    int Np = 100;

    // set the particle direction and position
    for (int n = 0; n < Np; ++n)
    {
        // sample a position WITHIN the geometry
        double x = 3.78 * p.rng().ran();
        double y = 3.78 * p.rng().ran();
        double z = 14.28 * p.rng().ran();

        double costheta = 1.0 - 2.0 * rng.ran();
        double phi      = profugus::constants::two_pi * rng.ran();
        double sintheta = sqrt(1.0 - costheta * costheta);

        double Ox = sintheta * cos(phi);
        double Oy = sintheta * sin(phi);
        double Oz = costheta;

        p.set_wt(1.0);

        // initialize geometry state
        geometry->initialize(Vector(x, y, z), Vector(Ox, Oy, Oz), geo);
        p.set_matid(geometry->matid(geo));
        EXPECT_TRUE(geometry->boundary_state(geo) ==
                    profugus::geometry::INSIDE);

        // initialize physics state
        physics->initialize(1.1, p);
        EXPECT_EQ(0, p.group());

        // transport the particle
        p.live();
        transporter.transport(p, bank);
        EXPECT_TRUE(!p.alive());

        // count up events
        if (p.event() == profugus::events::ESCAPE)
        {
            esc++;
        }
        else if (p.event() == profugus::events::ROULETTE_KILLED)
        {
            rk++;
        }
        else
        {
            ostringstream m;
            m << "Registered an impossible event, " << p.event();
            FAIL() << m.str();
        }
    }

    cout << "Events on external boundary-mesh:\n"
         << "\t" << esc << " escaping particles\n"
         << "\t" << rk  << " rouletted particles" << endl << endl;

    EXPECT_EQ(67, rk);
    EXPECT_EQ(33, esc);
}

//---------------------------------------------------------------------------//
//                 end of tstDomain_Transporter.cc
//---------------------------------------------------------------------------//
