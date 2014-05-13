//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstSource_Transporter.cc
 * \author Thomas M. Evans
 * \date   Tue May 13 13:54:23 2014
 * \brief  Source_Transporter unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Source_Transporter.hh"

#include "gtest/utils_gtest.hh"

#include <cmath>
#include <memory>

#include "comm/P_Stream.hh"
#include "comm/global.hh"
#include "utils/Definitions.hh"

#include "TransporterTestBase.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class DRSourceTransporterTest : public TransporterTestBase
{
    typedef TransporterTestBase Base;

  public:

    typedef profugus::Source_Transporter Transporter_t;

    void init_tallies()
    {
        // no tallies have been added
        tallier->build();
    }
};

//---------------------------------------------------------------------------//
// Replicated source

class DR_Source : public profugus::Source
{
    typedef profugus::Source Base;

  private:
    // number of particles per process
    size_type d_Np;
    size_type d_total;
    size_type d_running;

  public:
    DR_Source(SP_Geometry geometry, SP_Physics physics, SP_RNG_Control rcon)
        : Base(geometry, physics, rcon)
    {
        /* * */
    }

    SP_Particle get_particle()
    {
        Require (d_Np);

        using def::X; using def::Y; using def::Z;

        SP_Particle p(std::make_shared<Particle_t>());

        // get random number state for the particle
        auto rng = b_rng_control->rng();
        p->set_rng(rng);

        // sample position
        Space_Vector r;
        r[X] = 3.78  * p->rng().ran();
        r[Y] = 3.78  * p->rng().ran();
        r[Z] = 14.28 * p->rng().ran();

        // sample direction
        double costheta = 1.0 - 2.0 * rng.ran();
        double phi      = profugus::constants::two_pi * rng.ran();
        double sintheta = std::sqrt(1.0 - costheta * costheta);

        Space_Vector omega;
        omega[X] = sintheta * std::cos(phi);
        omega[Y] = sintheta * std::sin(phi);
        omega[Z] = costheta;

        // set position and direction
        p->set_wt(1.0);

        // initialize the physics and geometry state
        b_geometry->initialize(r, omega, p->geo_state());
        p->set_matid(b_geometry->matid(p->geo_state()));
        b_physics->initialize(1.1, *p);

        // make particle alive
        p->live();

        // subtract from total
        d_running--;

        return p;
    }

    bool empty() const
    {
        return !d_running;
    }

    size_type num_to_transport() const { return d_Np; }
    size_type num_to_transport_in_set() const { return d_Np; }
    size_type total_num_to_transport() const { return d_total; }

    void set_Np(size_type Np)
    {
        d_Np      = Np;
        d_running = d_Np;
        d_total   = d_Np * profugus::nodes();
    }

    size_type num_run() const { return d_Np - d_running; }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(DRSourceTransporterTest, source)
{
    // make the source
    std::shared_ptr<DR_Source> source(std::make_shared<DR_Source>(
                                          geometry, physics, rcon));
    source->set_Np(11);

    Source_t& base = *source;
    EXPECT_EQ(11, source->num_to_transport());
    EXPECT_EQ(nodes * 11, source->total_num_to_transport());
    EXPECT_EQ(0, source->num_run());

    int count = 0;
    while (!base.empty())
    {
        SP_Particle p = base.get_particle();
        EXPECT_TRUE(static_cast<bool>(p));
        count++;
    }

    EXPECT_EQ(11, count);
    EXPECT_EQ(11, source->num_run());
    EXPECT_EQ(11, source->num_to_transport());
    EXPECT_EQ(nodes * 11, source->total_num_to_transport());
}

//---------------------------------------------------------------------------//

TEST_F(DRSourceTransporterTest, Heuristic)
{
    db->set("mc_diag_frac", 0.2);

    // make the fixed source Transporter_t
    Transporter_t solver(db, geometry, physics);

    // set the variance reduction
    solver.set(var_red);

    // set the tally
    solver.set(tallier);

    // make the source
    std::shared_ptr<DR_Source> source(std::make_shared<DR_Source>(
                                          geometry, physics, rcon));
    source->set_Np(50);

    // assign the source
    solver.assign_source(source);

    // solve
    profugus::pcout << profugus::endl;
    solver.solve();
    profugus::pcout << profugus::endl;
}

//---------------------------------------------------------------------------//
//                 end of tstSource_Transporter.cc
//---------------------------------------------------------------------------//
