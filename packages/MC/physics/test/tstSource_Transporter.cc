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
#include "geometry/Geometry.hh"
#include "Shape.hh"
#include "Uniform_Source.hh"
#include "Box_Shape.hh"

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
// TESTS
//---------------------------------------------------------------------------//

TEST_F(DRSourceTransporterTest, source)
{
    // make a sampling shape (uniform)
    SP_Shape box(std::make_shared<profugus::Box_Shape>(
                     0.0, 2.52, 0.0, 2.52, 0.0, 14.28));

    // make the source
    int np = 11;
    db->set("Np", np);
    std::shared_ptr<Source_t> source(std::make_shared<Source_t>(
					 db,geometry, physics,box));

    Source_t& base = *source;
    EXPECT_EQ(11, source->num_to_transport());
    EXPECT_EQ(nodes * 11, source->total_num_to_transport());

    int count = 0;
    for ( int n = 0; n < np; ++n )
    {
        SP_Particle p = base.get_particle(n);
        EXPECT_TRUE(static_cast<bool>(p));
        count++;
    }

    EXPECT_EQ(11, count);
    EXPECT_EQ(11, source->num_to_transport());
    EXPECT_EQ(nodes * 11, source->total_num_to_transport());
}

//---------------------------------------------------------------------------//

TEST_F(DRSourceTransporterTest, Heuristic)
{
    db->set("mc_diag_frac", 0.2);
    db->set("Np", 11);

    // make the fixed source Transporter_t
    Transporter_t solver(db, geometry, physics);

    // set the variance reduction
    solver.set(var_red);

    // set the tally
    solver.set(tallier);

    // make a sampling shape (uniform)
    SP_Shape box(std::make_shared<profugus::Box_Shape>(
                     0.0, 2.52, 0.0, 2.52, 0.0, 14.28));
    // make the source
    std::shared_ptr<Source_t> source(std::make_shared<Source_t>(
					 db,geometry, physics,box));

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
