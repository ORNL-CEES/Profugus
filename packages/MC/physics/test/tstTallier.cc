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
#include "../Physics.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Tally Types
//---------------------------------------------------------------------------//

class A_Tally : public profugus::Pathlength_Tally
{
    typedef profugus::Pathlength_Tally Base;
  public:
    A_Tally(SP_Physics physics)
        : Base(physics, false)
    {
        set_name("a_pl_tally");
    }

    void accumulate(const double step, Particle_t &p) const { /* * */ }
};

//---------------------------------------------------------------------------//

class P_Tally : public profugus::Pathlength_Tally
{
    typedef profugus::Pathlength_Tally Base;
  public:
    P_Tally(SP_Physics physics)
        : Base(physics, false)
    {
        set_name("p_pl_tally");
    }

    void accumulate(const double step, Particle_t &p) const { /* * */ }
};

//---------------------------------------------------------------------------//

class Q_Tally : public profugus::Source_Tally
{
    typedef profugus::Source_Tally Base;
  public:
    Q_Tally(SP_Physics physics)
        : Base(physics, false)
    {
        set_name("q_src_tally");
    }

    void birth(const Particle_t &p) { /* * */ }
};

//---------------------------------------------------------------------------//

class S_Tally : public profugus::Source_Tally
{
    typedef profugus::Source_Tally Base;
  public:
    S_Tally(SP_Physics physics)
        : Base(physics, false)
    {
        set_name("s_src_tally");
    }

    void birth(const Particle_t &p) { /* * */ }
};

//---------------------------------------------------------------------------//

class C_Tally : public profugus::Compound_Tally
{
    typedef profugus::Compound_Tally Base;

  public:
    class C_PL_Tally : public profugus::Pathlength_Tally
    {
      public:
        C_PL_Tally(SP_Physics physics)
            : profugus::Pathlength_Tally(physics, false)
        {
            set_name("c_tally");
        }

        void accumulate(const double step, Particle_t &p) const { /* * */ }
    };

    class C_SRC_Tally : public profugus::Source_Tally
    {
      public:
        C_SRC_Tally(SP_Physics physics)
            : profugus::Source_Tally(physics, false)
        {
            set_name("c_tally");
        }

        void birth(const Particle_t &p) { /* * */ }
    };

  public:
    C_Tally(SP_Physics physics)
        : Base(physics, false)
    {
        set_name("c_tally");

        // make compound tallies
        b_pl_tally  = std::make_shared<C_PL_Tally>(physics);
        b_src_tally = std::make_shared<C_SRC_Tally>(physics);
    }
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class TallierTest : public testing::Test
{
  protected:
    typedef profugus::Tallier              Tallier_t;
    typedef Tallier_t::SP_Tally            SP_Tally;
    typedef profugus::Physics              Physics_t;
    typedef Physics_t::Geometry_t          Geometry_t;
    typedef Physics_t::Particle_t          Particle_t;
    typedef Physics_t::XS_t                XS_t;
    typedef Physics_t::RCP_XS              RCP_XS;
    typedef Physics_t::SP_Geometry         SP_Geometry;
    typedef Tallier_t::SP_Physics          SP_Physics;
    typedef Tallier_t::Pathlength_Tally_t  Pathlength_Tally_t;
    typedef Tallier_t::Source_Tally_t      Source_Tally_t;
    typedef Tallier_t::SP_Pathlength_Tally SP_Pathlength_Tally;
    typedef Tallier_t::SP_Source_Tally     SP_Source_Tally;

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
//---------------------------------------------------------------------------//
TEST_F(TallierTest, add)
{
    Tallier_t tallier;
    tallier.set(geometry, physics);

    // make tallies
    auto a(std::make_shared<A_Tally>(physics));
    auto p(std::make_shared<P_Tally>(physics));
    auto q(std::make_shared<Q_Tally>(physics));
    auto s(std::make_shared<S_Tally>(physics));
    auto c(std::make_shared<C_Tally>(physics));

    // add the tallies
    tallier.add_pathlength_tally(a);
    tallier.add_pathlength_tally(p);
    tallier.add_source_tally(q);
    tallier.add_source_tally(s);

    EXPECT_EQ(0, tallier.num_tallies());
    EXPECT_EQ(2, tallier.num_pathlength_tallies());
    EXPECT_EQ(2, tallier.num_source_tallies());
    EXPECT_EQ(0, tallier.num_compound_tallies());

    tallier.build();

    EXPECT_EQ(4, tallier.num_tallies());
    EXPECT_EQ(2, tallier.num_pathlength_tallies());
    EXPECT_EQ(2, tallier.num_source_tallies());
    EXPECT_EQ(0, tallier.num_compound_tallies());

    tallier.finalize(1);
    tallier.reset();

    EXPECT_EQ(0, tallier.num_tallies());
    EXPECT_EQ(2, tallier.num_pathlength_tallies());
    EXPECT_EQ(2, tallier.num_source_tallies());
    EXPECT_EQ(0, tallier.num_compound_tallies());

    // add duplicates
    tallier.add_pathlength_tally(a);
    tallier.add_source_tally(q);

    EXPECT_EQ(0, tallier.num_tallies());
    EXPECT_EQ(3, tallier.num_pathlength_tallies());
    EXPECT_EQ(3, tallier.num_source_tallies());
    EXPECT_EQ(0, tallier.num_compound_tallies());

    tallier.build();

    EXPECT_EQ(4, tallier.num_tallies());
    EXPECT_EQ(2, tallier.num_pathlength_tallies());
    EXPECT_EQ(2, tallier.num_source_tallies());
    EXPECT_EQ(0, tallier.num_compound_tallies());

    tallier.finalize(1);
    tallier.reset();

    EXPECT_EQ(0, tallier.num_tallies());
    EXPECT_EQ(2, tallier.num_pathlength_tallies());
    EXPECT_EQ(2, tallier.num_source_tallies());
    EXPECT_EQ(0, tallier.num_compound_tallies());

    // add compound tallies
    tallier.add_compound_tally(c);

    EXPECT_EQ(0, tallier.num_tallies());
    EXPECT_EQ(3, tallier.num_pathlength_tallies());
    EXPECT_EQ(3, tallier.num_source_tallies());
    EXPECT_EQ(1, tallier.num_compound_tallies());

    tallier.build();

    EXPECT_EQ(7, tallier.num_tallies());
    EXPECT_EQ(3, tallier.num_pathlength_tallies());
    EXPECT_EQ(3, tallier.num_source_tallies());
    EXPECT_EQ(1, tallier.num_compound_tallies());
}

//---------------------------------------------------------------------------//
//                 end of tstTallier.cc
//---------------------------------------------------------------------------//
