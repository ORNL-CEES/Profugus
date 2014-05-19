//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstTallier.cc
 * \author Thomas M. Evans
 * \date   Fri May 16 13:50:28 2014
 * \brief  Tallier test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Tallier.hh"
#include "../Keff_Tally.hh"
#include "../Physics.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class TallierTest : public testing::Test
{
  protected:
    typedef profugus::Tallier      Tallier_t;
    typedef Tallier_t::SP_Tally    SP_Tally;
    typedef profugus::Physics      Physics_t;
    typedef Physics_t::Geometry_t  Geometry_t;
    typedef profugus::Keff_Tally   Keff_Tally_t;
    typedef Physics_t::Particle_t  Particle_t;
    typedef Physics_t::XS_t        XS_t;
    typedef Physics_t::RCP_XS      RCP_XS;
    typedef Physics_t::SP_Geometry SP_Geometry;
    typedef Tallier_t::SP_Physics  SP_Physics;

    typedef Teuchos::ParameterList        ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t> RCP_Std_DB;

  protected:
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();

        db = Teuchos::rcp(new ParameterList_t("test"));

        build_geometry();
        build_physics();
    }

    void build_geometry()
    {
        typedef Geometry_t::SP_Array SP_Core;
        typedef Geometry_t::Array_t  Core_t;
        typedef Core_t::SP_Object    SP_Lattice;
        typedef Core_t::Object_t     Lattice_t;
        typedef Lattice_t::SP_Object SP_Pin_Cell;
        typedef Lattice_t::Object_t  Pin_Cell_t;

        // make boxes (10x10x20)
        SP_Pin_Cell p0(std::make_shared<Pin_Cell_t>(0, 10.0, 20.0));
        SP_Pin_Cell p1(std::make_shared<Pin_Cell_t>(1, 10.0, 20.0));

        // make lattice
        /*
          1 0
          0 1
         */
        SP_Lattice lat(std::make_shared<Lattice_t>(2, 2, 1, 2));

        // assign pins
        lat->assign_object(p0, 0); // mat 0 box
        lat->assign_object(p1, 1); // mat 1 box

        // arrange pin-cells in lattice
        lat->id(0, 0, 0) = 0; // mat 0 box
        lat->id(1, 0, 0) = 1; // mat 1 box
        lat->id(0, 1, 0) = 1; // mat 1 box
        lat->id(1, 1, 0) = 0; // mat 0 box

        // complete lattice
        lat->complete(0.0, 0.0, 0.0);

        // make core
        SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));
        core->assign_object(lat, 0);
        core->complete(0.0, 0.0, 0.0);

        geometry = std::make_shared<Geometry_t>(core);

    }

    void build_physics()
    {
        const int ng = 3;

        RCP_XS xs(Teuchos::rcp(new XS_t()));
        xs->set(0, ng);

        // make group boundaries
        XS_t::OneDArray nbnd(4, 0.0);
        nbnd[0] = 100.0; nbnd[1] = 1.0; nbnd[2] = 0.01; nbnd[3] = 0.0001;
        xs->set_bounds(nbnd);

        double t1[] = {1.1, 1.6, 2.9};
        double t2[] = {10.0, 11.3, 16.2};

        XS_t::OneDArray tot1(std::begin(t1), std::end(t1));
        XS_t::OneDArray tot2(std::begin(t2), std::end(t2));

        xs->add(0, XS_t::TOTAL, tot1);
        xs->add(1, XS_t::TOTAL, tot2);

        double s1[][3] = {{0.7, 0.0, 0.0},
                          {0.2, 0.3, 0.0},
                          {0.1, 0.7, 1.9}};

        double s2[][3] = {{2.7, 0.0, 0.0},
                          {2.2, 2.3, 0.0},
                          {2.1, 2.7, 3.9}};

        XS_t::TwoDArray sct1(ng, ng, 0.0);
        XS_t::TwoDArray sct2(ng, ng, 0.0);

        for (int g = 0; g < 3; ++g)
        {
            for (int gp = 0; gp < 3; ++gp)
            {
                sct1(g, gp) = s1[g][gp];
                sct2(g, gp) = s2[g][gp];
            }
        }

        xs->add(0, 0, sct1);
        xs->add(1, 0, sct2);

        double c2[] = {0.4, 0.6, 0.0};
        double f2[] = {3.2, 4.2, 0.0};
        double n2[] = {2.4*3.2, 2.4*4.2, 0.0};

        XS_t::OneDArray chi2(std::begin(c2), std::end(c2));
        XS_t::OneDArray fis2(std::begin(f2), std::end(f2));
        XS_t::OneDArray nuf2(std::begin(n2), std::end(n2));

        xs->add(1, XS_t::CHI, chi2);
        xs->add(1, XS_t::SIG_F, fis2);
        xs->add(1, XS_t::NU_SIG_F, nuf2);

        xs->complete();

        physics = std::make_shared<Physics_t>(db, xs);
        physics->set_geometry(geometry);
    }

  protected:
    RCP_Std_DB  db;
    SP_Geometry geometry;
    SP_Physics  physics;

    int node;
    int nodes;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(TallierTest, keff_tally)
{
    Tallier_t tallier;
    tallier.set(geometry, physics);

    EXPECT_EQ(0, tallier.num_tallies());
    EXPECT_EQ(0, tallier.num_pathlength_tallies());
    EXPECT_EQ(0, tallier.num_source_tallies());
    EXPECT_FALSE(tallier.is_built());
    EXPECT_FALSE(tallier.is_finalized());

    // make a keff tally
    auto keff(std::make_shared<Keff_Tally_t>(1.0, physics));

    // add tallies
    tallier.add_pathlength_tally(keff);
    EXPECT_EQ(0, tallier.num_tallies());
    EXPECT_EQ(1, tallier.num_pathlength_tallies());
    EXPECT_EQ(0, tallier.num_source_tallies());
    EXPECT_FALSE(tallier.is_built());
    EXPECT_FALSE(tallier.is_finalized());

    // build the tallier
    tallier.build();
    EXPECT_EQ(1, tallier.num_tallies());
    EXPECT_EQ(1, tallier.num_pathlength_tallies());
    EXPECT_EQ(0, tallier.num_source_tallies());
    EXPECT_TRUE(tallier.is_built());
    EXPECT_FALSE(tallier.is_finalized());

    // Make a particle
    auto p = std::make_shared<Particle_t>();

    // Create reference variables
    double ref_k     = 0.0;
    double ref_kavg  = 0.0;
    double ref_kmom1 = 0.0;
    double ref_kmom2 = 0.0;
    double ref_kvar  = 0.0;

    /*** Begin INACTIVE CYCLE 1 ***/
    tallier.begin_cycle();
    EXPECT_EQ(0.0, keff->latest());

    // Tally a particle
    p->set_wt(0.65);
    p->set_matid(1);
    p->set_group(1);
    tallier.path_length(1.4, *p);
    ref_k = 1.4 * 0.65 * 2.4 * 4.2;
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in group 2
    p->set_group(2);  // no fission in group 2
    tallier.path_length(1.8, *p);
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in group 0
    p->set_group(0);
    p->set_wt(0.4);
    tallier.path_length(0.6, *p);
    ref_k += 0.4 * 0.6 * 2.4 * 3.2;
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in material 0
    p->set_matid(0);
    p->set_wt(0.55);
    tallier.path_length(0.2, *p);
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    tallier.end_history();

    // End the cycle
    tallier.end_cycle(3.0);
    ref_k = nodes * ref_k / 3.0;
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    /*** Begin INACTIVE CYCLE 2 ***/
    tallier.begin_cycle();
    EXPECT_EQ(0.0, keff->latest());

    // Tally more particles
    p->set_wt(0.55);
    p->set_matid(1);
    p->set_group(1);
    tallier.path_length(1.4, *p);
    ref_k = 1.4 * 0.55 * 2.4 * 4.2;
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in group 2
    p->set_group(2);  // no fission in group 2
    tallier.path_length(1.8, *p);
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in group 0
    p->set_group(0);
    p->set_wt(0.2);
    tallier.path_length(0.6, *p);
    ref_k += 0.2 * 0.6 * 2.4 * 3.2;
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in material 0
    p->set_matid(0);
    p->set_wt(0.55);
    tallier.path_length(0.2, *p);
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    tallier.end_history();

    // End the cycle
    tallier.end_cycle(3.0);
    ref_k = nodes * ref_k / 3.0;
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    /*** Begin ACTIVE CYCLE 1 ***/

    tallier.begin_active_cycles();
    tallier.begin_cycle();
    EXPECT_EQ(0.0, keff->latest());
    EXPECT_EQ(0.0, keff->keff_sum());
    EXPECT_EQ(0.0, keff->keff_sum_sq());

    // Tally some particles
    p->set_wt(0.65);
    p->set_matid(1);
    p->set_group(1);
    tallier.path_length(1.4, *p);
    ref_k = 1.4 * 0.65 * 2.4 * 4.2;
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in group 2
    p->set_group(2);  // no fission in group 2
    tallier.path_length(1.8, *p);
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in group 0
    p->set_group(0);
    p->set_wt(0.4);
    tallier.path_length(0.6, *p);
    ref_k += 0.4 * 0.6 * 2.4 * 3.2;
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in material 0
    p->set_matid(0);
    p->set_wt(0.55);
    tallier.path_length(0.2, *p);
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    tallier.end_history();

    // End the cycle
    tallier.end_cycle(3.0);

    ref_k = nodes * ref_k / 3.0;
    ref_kmom1 += ref_k;
    ref_kmom2 += ref_k * ref_k;
    ref_kavg = ref_kmom1 / 1.0;

    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);
    EXPECT_SOFTEQ(ref_kavg, keff->mean(), 1.0e-6);

    /*** Begin ACTIVE CYCLE 2 ***/

    tallier.begin_cycle();
    EXPECT_EQ(0.0, keff->latest());
    EXPECT_NE(0.0, keff->keff_sum());
    EXPECT_NE(0.0, keff->keff_sum_sq());

    // Tally some more particles
    p->set_wt(0.55);
    p->set_matid(1);
    p->set_group(1);
    tallier.path_length(1.4, *p);
    ref_k = 1.4 * 0.55 * 2.4 * 4.2;
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in group 2
    p->set_group(2);  // no fission in group 2
        tallier.path_length(1.8, *p);
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in group 0
    p->set_group(0);
    p->set_wt(0.2);
    tallier.path_length(0.6, *p);
    ref_k += 0.2 * 0.6 * 2.4 * 3.2;
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in material 0
    p->set_matid(0);
    p->set_wt(0.55);
    tallier.path_length(0.2, *p);
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    tallier.end_history();

    // End the cycle
    tallier.end_cycle(3.0);

    ref_k = nodes * ref_k / 3.0;
    ref_kmom1 += ref_k;
    ref_kmom2 += ref_k * ref_k;
    ref_kavg = ref_kmom1 / 2.0;
    ref_kvar = (ref_kmom2 - 2.0 * ref_kavg * ref_kavg) / 2.0 / 1.0;

    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);
    EXPECT_SOFTEQ(ref_kavg, keff->mean(), 1.0e-6);
    EXPECT_SOFTEQ(ref_kvar, keff->variance(), 1.0e-6);

    /*** Begin ACTIVE CYCLE 3 ***/

    tallier.begin_cycle();
    EXPECT_EQ(0.0, keff->latest());
    EXPECT_NE(0.0, keff->keff_sum());
    EXPECT_NE(0.0, keff->keff_sum_sq());

    // Tally some more particles
    p->set_wt(0.8);
    p->set_matid(1);
    p->set_group(1);
    tallier.path_length(1.4, *p);
    ref_k = 1.4 * 0.8 * 2.4 * 4.2;
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in group 2
    p->set_group(2);  // no fission in group 2
    tallier.path_length(1.8, *p);
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in group 0
    p->set_group(0);
    p->set_wt(0.65);
    tallier.path_length(0.6, *p);
    ref_k += 0.65 * 0.6 * 2.4 * 3.2;
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    // Put particle in material 0
    p->set_matid(0);
    p->set_wt(0.55);
    tallier.path_length(0.2, *p);
    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);

    tallier.end_history();

    // End the cycle
    tallier.end_cycle(3.0);

    EXPECT_FALSE(tallier.is_finalized());
    tallier.finalize(15.0);
    EXPECT_TRUE(tallier.is_finalized());

    ref_k = nodes * ref_k / 3.0;
    ref_kmom1 += ref_k;
    ref_kmom2 += ref_k * ref_k;
    ref_kavg = ref_kmom1 / 3.0;
    ref_kvar = (ref_kmom2 - 3.0 * ref_kavg * ref_kavg) / 3.0 / 2.0;

    EXPECT_SOFTEQ(ref_k, keff->latest(), 1.0e-6);
    EXPECT_SOFTEQ(ref_kavg, keff->mean(), 1.0e-6);
    EXPECT_SOFTEQ(ref_kvar, keff->variance(), 1.0e-6);

    tallier.reset();

    EXPECT_EQ(0, tallier.num_tallies());
    EXPECT_EQ(1, tallier.num_pathlength_tallies());
    EXPECT_EQ(0, tallier.num_source_tallies());
    EXPECT_FALSE(tallier.is_built());
    EXPECT_FALSE(tallier.is_finalized());

    EXPECT_EQ(0., keff->keff_sum());
    EXPECT_EQ(0., keff->keff_sum_sq());
    EXPECT_EQ(0, keff->cycle_count());
}

//---------------------------------------------------------------------------//

TEST_F(TallierTest, swap)
{
    Tallier_t tallier, inactive_tallier;
    tallier.set(geometry, physics);
    inactive_tallier.set(geometry, physics);

    // make a keff tally
    auto keff(std::make_shared<Keff_Tally_t>(1.0, physics));

    // add tallies
    tallier.add_pathlength_tally(keff);

    // build the talliers
    tallier.build();
    inactive_tallier.build();

    EXPECT_EQ(1, tallier.num_tallies());
    EXPECT_EQ(1, tallier.num_pathlength_tallies());
    EXPECT_EQ(0, tallier.num_source_tallies());
    EXPECT_TRUE(tallier.is_built());
    EXPECT_FALSE(tallier.is_finalized());

    EXPECT_EQ(0, inactive_tallier.num_tallies());
    EXPECT_EQ(0, inactive_tallier.num_pathlength_tallies());
    EXPECT_EQ(0, inactive_tallier.num_source_tallies());
    EXPECT_TRUE(inactive_tallier.is_built());
    EXPECT_FALSE(inactive_tallier.is_finalized());

    // swap the talliers
    swap(tallier, inactive_tallier);

    EXPECT_EQ(1, inactive_tallier.num_tallies());
    EXPECT_EQ(1, inactive_tallier.num_pathlength_tallies());
    EXPECT_EQ(0, inactive_tallier.num_source_tallies());
    EXPECT_TRUE(inactive_tallier.is_built());
    EXPECT_FALSE(inactive_tallier.is_finalized());

    EXPECT_EQ(0, tallier.num_tallies());
    EXPECT_EQ(0, tallier.num_pathlength_tallies());
    EXPECT_EQ(0, tallier.num_source_tallies());
    EXPECT_TRUE(tallier.is_built());
    EXPECT_FALSE(tallier.is_finalized());

    // swap the talliers
    swap(tallier, inactive_tallier);

    EXPECT_EQ(1, tallier.num_tallies());
    EXPECT_EQ(1, tallier.num_pathlength_tallies());
    EXPECT_EQ(0, tallier.num_source_tallies());
    EXPECT_TRUE(tallier.is_built());
    EXPECT_FALSE(tallier.is_finalized());

    EXPECT_EQ(0, inactive_tallier.num_tallies());
    EXPECT_EQ(0, inactive_tallier.num_pathlength_tallies());
    EXPECT_EQ(0, inactive_tallier.num_source_tallies());
    EXPECT_TRUE(inactive_tallier.is_built());
    EXPECT_FALSE(inactive_tallier.is_finalized());
}

//---------------------------------------------------------------------------//
//                 end of tstTallier.cc
//---------------------------------------------------------------------------//
