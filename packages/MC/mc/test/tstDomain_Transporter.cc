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
    typedef profugus::Domain_Transporter       Transporter;
    typedef Physics_t::Fission_Site_Container  Fission_Site_Container;
    typedef shared_ptr<Fission_Site_Container> SP_Fission_Sites;

  protected:

    void init_vr()
    {
        db->set("weight_cutoff", 0.001);
        var_red = std::make_shared<profugus::VR_Roulette>(db);
    }

    void init_tallies()
    {
        // no tallies have been added
        tallier->build();
    }
};

//---------------------------------------------------------------------------//

class Reflecting_Domain_TransporterTest : public TransporterTestBase
{
  protected:
    typedef profugus::Domain_Transporter Transporter;

  protected:

    void init_vr()
    {
        db->set("weight_cutoff", 0.001);
        var_red = std::make_shared<profugus::VR_Roulette>(db);
    }

    void init_tallies()
    {
        // no tallies have been added
        tallier->build();
    }

    void init_geometry()
    {
        typedef Geometry_t::SP_Array SP_Core;
        typedef Geometry_t::Array_t  Core_t;
        typedef Core_t::SP_Object    SP_Lattice;
        typedef Core_t::Object_t     Lattice_t;
        typedef Lattice_t::SP_Object SP_Pin_Cell;
        typedef Lattice_t::Object_t  Pin_Cell_t;
        typedef Core_t::Vec_Int      Vec_Int;

        // make pin cells
        SP_Pin_Cell p1(std::make_shared<Pin_Cell_t>(1, 0.54, 3, 1.26, 14.28));
        SP_Pin_Cell p2(std::make_shared<Pin_Cell_t>(2, 0.54, 3, 1.26, 14.28));

        // make lattice
        SP_Lattice lat(std::make_shared<Lattice_t>(3, 3, 1, 2));

        // assign pins
        lat->assign_object(p1, 0); // fuel pins
        lat->assign_object(p2, 1); // guide tube

        // arrange pin-cells in lattice
        lat->id(0, 0, 0) = 0; // fuel pin
        lat->id(1, 0, 0) = 0; // fuel pin
        lat->id(2, 0, 0) = 0; // fuel pin
        lat->id(0, 1, 0) = 0; // fuel pin
        lat->id(1, 1, 0) = 1; // guide tube
        lat->id(2, 1, 0) = 0; // fuel pin
        lat->id(0, 2, 0) = 0; // fuel pin
        lat->id(1, 2, 0) = 0; // fuel pin
        lat->id(2, 2, 0) = 0; // fuel pin

        // complete lattice
        lat->complete(0.0, 0.0, 0.0);

        Check(profugus::soft_equiv(lat->pitch(def::X), 3.78));
        Check(profugus::soft_equiv(lat->pitch(def::Y), 3.78));
        Check(profugus::soft_equiv(lat->height(), 14.28));

        // make core
        SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));
        core->assign_object(lat, 0);
        core->set_reflecting(Vec_Int(6, 1));
        core->complete(0.0, 0.0, 0.0);

        geometry = std::make_shared<Geometry_t>(core);
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

TEST_F(Reflecting_Domain_TransporterTest, transport)
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

    EXPECT_EQ(100, rk);
    EXPECT_EQ(0, esc);
}

//---------------------------------------------------------------------------//

TEST_F(Domain_TransporterTest, fission_sampling)
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

    // make a fission site container
    Transporter::SP_Fission_Sites fis_sites(
        std::make_shared<Fission_Site_Container>());

    // set the fission site container in the transport so that fission sites
    // get sampled
    transporter.set(fis_sites, 1.0);

    // set the particle direction and position
    for (int n = 0; n < Np; ++n)
    {
        // sample a position WITHIN the geometry mesh
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
        EXPECT_TRUE(geometry->boundary_state(geo)
                    == profugus::geometry::INSIDE);

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

    // make sure all the sites are in fission
    int nfp = 0;
    Particle_t::Geo_State_t &gs = p.geo_state();
    for (int n = 0; n < fis_sites->size(); ++n)
    {
        auto &sites = *fis_sites;
        Vector r    = physics->fission_site(sites[n]);
        nfp        += 1;

        // make sure each site is in a fissionable material
        geometry->initialize(r, Vector(1.0, 0.0, 0.0), gs);
        EXPECT_EQ(1, geometry->matid(gs));
        p.set_matid(geometry->matid(gs));

        // check group, etc....
        physics->initialize_fission(sites[n], p);
        EXPECT_EQ(4.0, physics->total(profugus::physics::FISSION, p));
    }

    // Number of fission samples
    cout << "Sampled " << fis_sites->size() << " separate fission sites "
         << "containing " << nfp << " Fission source particles" << endl << endl;

    EXPECT_EQ(74, rk);
    EXPECT_EQ(26, esc);
}

//---------------------------------------------------------------------------//
//                 end of tstDomain_Transporter.cc
//---------------------------------------------------------------------------//
