//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/test/tstXS_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Feb 05 20:07:26 2014
 * \brief  XS_Builder unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../XS_Builder.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class XS_Builder_Test : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture
    typedef profugus::XS_Builder   XS_Builder;
    typedef XS_Builder::XS_t       XS_t;
    typedef XS_Builder::RCP_XS     RCP_XS;
    typedef XS_Builder::Vec_Str    Vec_Str;
    typedef XS_Builder::Matid_Map  Matid_Map;
    typedef XS_Builder::std_string std_string;
    typedef XS_t::Vec_Int          Vec_Int;
    typedef XS_t::Vector           Vector;
    typedef XS_t::Matrix           Matrix;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();
    }

  protected:
    // >>> Data that get re-initialized between tests

    int node, nodes;
    XS_Builder builder;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(XS_Builder_Test, load_3GP0)
{
    builder.open_and_broadcast("xs3GP0.xml");
    EXPECT_EQ(3, builder.num_groups());
    EXPECT_EQ(0, builder.pn_order());

    const Vec_Str &mats = builder.materials();
    EXPECT_EQ(2, mats.size());
    EXPECT_EQ("mat 0", mats[0]);
    EXPECT_EQ("mat 10", mats[1]);

    // read all of the materials
    Matid_Map map;
    map.insert(Matid_Map::value_type(1, std_string("mat 0")));
    map.insert(Matid_Map::value_type(2, std_string("mat 10")));
    map.complete();

    // build the cross sections
    builder.build(map);
    RCP_XS xs = builder.get_xs();
    EXPECT_FALSE(xs.is_null());

    // test the xs
    EXPECT_EQ(0, xs->pn_order());
    EXPECT_EQ(3, xs->num_groups());

    Vec_Int matids;
    xs->get_matids(matids);
    EXPECT_EQ(2, matids.size());
    EXPECT_EQ(1, matids[0]);
    EXPECT_EQ(2, matids[1]);

    // test bounds
    const auto &bounds = xs->bounds();
    EXPECT_EQ(4, bounds.length());
    EXPECT_EQ(0.0, bounds(0));
    EXPECT_EQ(0.0, bounds(1));
    EXPECT_EQ(0.0, bounds(2));
    EXPECT_EQ(0.0, bounds(3));

    // test velocities
    const auto &vel = xs->velocities();
    EXPECT_EQ(3, vel.length());
    EXPECT_EQ(0.0, vel(0));
    EXPECT_EQ(0.0, vel(1));
    EXPECT_EQ(0.0, vel(2));

    // material 1
    {
        const Vector &sigt = xs->vector(1, XS_t::TOTAL);
        EXPECT_EQ(3, sigt.length());

        EXPECT_EQ(1.0, sigt(0));
        EXPECT_EQ(2.0, sigt(1));
        EXPECT_EQ(3.0, sigt(2));

        for (int t = 1; t < XS_t::END_XS_TYPES; ++t)
        {
            const Vector &sig = xs->vector(1, t);
            EXPECT_EQ(3, sig.length());
            for (int g = 0; g < 3; ++g)
            {
                EXPECT_EQ(0.0, sig(g));
            }
        }

        const Matrix &p0 = xs->matrix(1, 0);
        EXPECT_EQ(3, p0.numRows());
        EXPECT_EQ(3, p0.numCols());

        EXPECT_EQ(0.6, p0(0, 0));
        EXPECT_EQ(0.0, p0(0, 1));
        EXPECT_EQ(0.0, p0(0, 2));

        EXPECT_EQ(0.2, p0(1, 0));
        EXPECT_EQ(1.8, p0(1, 1));
        EXPECT_EQ(0.1, p0(1, 2));

        EXPECT_EQ(0.1, p0(2, 0));
        EXPECT_EQ(0.1, p0(2, 1));
        EXPECT_EQ(2.5, p0(2, 2));
    }

    // material 2
    {
        const Vector &sigt = xs->vector(2, XS_t::TOTAL);
        EXPECT_EQ(3, sigt.length());

        EXPECT_EQ(10.0, sigt(0));
        EXPECT_EQ(11.0, sigt(1));
        EXPECT_EQ(12.0, sigt(2));

        const Vector &sigf = xs->vector(2, XS_t::SIG_F);
        EXPECT_EQ(3, sigf.length());

        EXPECT_EQ(3.0, sigf(0));
        EXPECT_EQ(4.5, sigf(1));
        EXPECT_EQ(5.0, sigf(2));

        const Vector &nusigf = xs->vector(2, XS_t::NU_SIG_F);
        EXPECT_EQ(3, nusigf.length());

        EXPECT_SOFTEQ(2.4*3.0, nusigf(0), 1.e-12);
        EXPECT_SOFTEQ(2.4*4.5, nusigf(1), 1.e-12);
        EXPECT_SOFTEQ(2.4*5.0, nusigf(2), 1.e-12);

        const Vector &chi = xs->vector(2, XS_t::CHI);
        EXPECT_EQ(3, chi.length());

        EXPECT_EQ(0.9, chi(0));
        EXPECT_EQ(0.1, chi(1));
        EXPECT_EQ(0.0, chi(2));

        const Matrix &p0 = xs->matrix(2, 0);
        EXPECT_EQ(3, p0.numRows());
        EXPECT_EQ(3, p0.numCols());

        EXPECT_EQ(1.7, p0(0, 0));
        EXPECT_EQ(0.0, p0(0, 1));
        EXPECT_EQ(0.0, p0(0, 2));

        EXPECT_EQ(0.4, p0(1, 0));
        EXPECT_EQ(5.1, p0(1, 1));
        EXPECT_EQ(0.4, p0(1, 2));

        EXPECT_EQ(0.3, p0(2, 0));
        EXPECT_EQ(1.2, p0(2, 1));
        EXPECT_EQ(3.0, p0(2, 2));
    }
}

//---------------------------------------------------------------------------//

TEST_F(XS_Builder_Test, load_5GP1)
{
    builder.open_and_broadcast("xs5GP1.xml");
    EXPECT_EQ(5, builder.num_groups());
    EXPECT_EQ(1, builder.pn_order());

    const Vec_Str &mats = builder.materials();
    EXPECT_EQ(3, mats.size());
    EXPECT_EQ("mat 1", mats[0]);
    EXPECT_EQ("mat 4", mats[1]);
    EXPECT_EQ("mat 5", mats[2]);

    // read 2 materials
    Matid_Map map;
    map.insert(Matid_Map::value_type(1, std_string("mat 4")));
    map.insert(Matid_Map::value_type(2, std_string("mat 5")));
    map.complete();

    // build the cross sections
    builder.build(map);
    RCP_XS xs = builder.get_xs();
    EXPECT_FALSE(xs.is_null());

    // test the xs
    EXPECT_EQ(1, xs->pn_order());
    EXPECT_EQ(5, xs->num_groups());
    EXPECT_EQ(2, xs->num_mat());

    EXPECT_TRUE(xs->has(1));
    EXPECT_TRUE(xs->has(2));

    // test bounds
    const auto &bounds = xs->bounds();
    EXPECT_EQ(6, bounds.length());
    EXPECT_EQ(2.0e7,  bounds(0));
    EXPECT_EQ(1.0e6,  bounds(1));
    EXPECT_EQ(1.0e5,  bounds(2));
    EXPECT_EQ(1.0e4,  bounds(3));
    EXPECT_EQ(1.0e3,  bounds(4));
    EXPECT_EQ(1.0e-5, bounds(5));

    // material 2
    {
        const Vector &sigt = xs->vector(2, XS_t::TOTAL);
        EXPECT_EQ(5, sigt.length());

        EXPECT_EQ(5, sigt(0));
        EXPECT_EQ(7.5, sigt(1));
        EXPECT_EQ(8, sigt(2));
        EXPECT_EQ(10, sigt(3));
        EXPECT_EQ(26.5, sigt(4));

        const Vector &sigf = xs->vector(2, XS_t::SIG_F);
        EXPECT_EQ(5, sigf.length());

        EXPECT_EQ(1.5, sigf(0));
        EXPECT_EQ(2.25, sigf(1));
        EXPECT_EQ(2.4, sigf(2));
        EXPECT_EQ(3.0, sigf(3));
        EXPECT_EQ(7.95, sigf(4));

        const Vector &nusigf = xs->vector(2, XS_t::NU_SIG_F);
        EXPECT_EQ(5, nusigf.length());

        EXPECT_EQ(3.6, nusigf(0));
        EXPECT_EQ(5.4, nusigf(1));
        EXPECT_EQ(5.76, nusigf(2));
        EXPECT_EQ(7.2, nusigf(3));
        EXPECT_EQ(19.08, nusigf(4));

        const Vector &chi = xs->vector(2, XS_t::CHI);
        EXPECT_EQ(5, chi.length());

        EXPECT_EQ(0.6, chi(0));
        EXPECT_EQ(0.3, chi(1));
        EXPECT_EQ(0.1, chi(2));
        EXPECT_EQ(0.0, chi(3));
        EXPECT_EQ(0.0, chi(4));

        const Matrix &p0 = xs->matrix(2, 0);
        const Matrix &p1 = xs->matrix(2, 1);
        EXPECT_EQ(5, p0.numRows());
        EXPECT_EQ(5, p0.numCols());
        EXPECT_EQ(5, p1.numRows());
        EXPECT_EQ(5, p1.numCols());

        EXPECT_EQ(0.8,  p0(0, 0));
        EXPECT_EQ(0.0,  p0(0, 1));
        EXPECT_EQ(0.0,  p0(0, 2));
        EXPECT_EQ(0.0,  p0(0, 3));
        EXPECT_EQ(0.0,  p0(0, 4));

        EXPECT_EQ(0.6,  p0(1, 0));
        EXPECT_EQ(1.5,  p0(1, 1));
        EXPECT_EQ(0.0,  p0(1, 2));
        EXPECT_EQ(0.0,  p0(1, 3));
        EXPECT_EQ(0.0,  p0(1, 4));

        EXPECT_EQ(0.2,  p0(2, 0));
        EXPECT_EQ(0.9,  p0(2, 1));
        EXPECT_EQ(1.92, p0(2, 2));
        EXPECT_EQ(0.0,  p0(2, 3));
        EXPECT_EQ(0.0,  p0(2, 4));

        EXPECT_EQ(0.2,  p0(3, 0));
        EXPECT_EQ(0.3,  p0(3, 1));
        EXPECT_EQ(0.96, p0(3, 2));
        EXPECT_EQ(2.8,  p0(3, 3));
        EXPECT_EQ(1.06, p0(3, 4));

        EXPECT_EQ(0.2,  p0(4, 0));
        EXPECT_EQ(0.3,  p0(4, 1));
        EXPECT_EQ(0.32, p0(4, 2));
        EXPECT_EQ(1.2,  p0(4, 3));
        EXPECT_EQ(9.54, p0(4, 4));

        for (int g = 0; g < 5; ++g)
        {
            for (int gp = 0; gp < 5; ++gp)
            {
                EXPECT_SOFTEQ(p0(g, gp) * 0.01, p1(g ,gp), 1.e-12);
            }
        }
    }
}

//---------------------------------------------------------------------------//

TEST_F(XS_Builder_Test, load_5GP1_P0)
{
    builder.open_and_broadcast("xs5GP1.xml");
    EXPECT_EQ(5, builder.num_groups());
    EXPECT_EQ(1, builder.pn_order());

    const Vec_Str &mats = builder.materials();
    EXPECT_EQ(3, mats.size());
    EXPECT_EQ("mat 1", mats[0]);
    EXPECT_EQ("mat 4", mats[1]);
    EXPECT_EQ("mat 5", mats[2]);

    // read 2 materials
    Matid_Map map;
    map.insert(Matid_Map::value_type(1, std_string("mat 4")));
    map.insert(Matid_Map::value_type(2, std_string("mat 5")));
    map.complete();

    // build the cross sections
    builder.build(map, 0, 0, 4);
    RCP_XS xs = builder.get_xs();
    EXPECT_FALSE(xs.is_null());

    // test the xs
    EXPECT_EQ(0, xs->pn_order());
    EXPECT_EQ(5, xs->num_groups());
    EXPECT_EQ(2, xs->num_mat());

    EXPECT_TRUE(xs->has(1));
    EXPECT_TRUE(xs->has(2));

    // material 2
    {
        const Vector &sigt = xs->vector(2, XS_t::TOTAL);
        EXPECT_EQ(5, sigt.length());

        EXPECT_EQ(5, sigt(0));
        EXPECT_EQ(7.5, sigt(1));
        EXPECT_EQ(8, sigt(2));
        EXPECT_EQ(10, sigt(3));
        EXPECT_EQ(26.5, sigt(4));

        const Vector &sigf = xs->vector(2, XS_t::SIG_F);
        EXPECT_EQ(5, sigf.length());

        EXPECT_EQ(1.5, sigf(0));
        EXPECT_EQ(2.25, sigf(1));
        EXPECT_EQ(2.4, sigf(2));
        EXPECT_EQ(3.0, sigf(3));
        EXPECT_EQ(7.95, sigf(4));

        const Vector &nusigf = xs->vector(2, XS_t::NU_SIG_F);
        EXPECT_EQ(5, nusigf.length());

        EXPECT_EQ(3.6, nusigf(0));
        EXPECT_EQ(5.4, nusigf(1));
        EXPECT_EQ(5.76, nusigf(2));
        EXPECT_EQ(7.2, nusigf(3));
        EXPECT_EQ(19.08, nusigf(4));

        const Vector &chi = xs->vector(2, XS_t::CHI);
        EXPECT_EQ(5, chi.length());

        EXPECT_EQ(0.6, chi(0));
        EXPECT_EQ(0.3, chi(1));
        EXPECT_EQ(0.1, chi(2));
        EXPECT_EQ(0.0, chi(3));
        EXPECT_EQ(0.0, chi(4));

        const Matrix &p0 = xs->matrix(2, 0);
        EXPECT_EQ(5, p0.numRows());
        EXPECT_EQ(5, p0.numCols());

        EXPECT_EQ(0.8,  p0(0, 0));
        EXPECT_EQ(0.0,  p0(0, 1));
        EXPECT_EQ(0.0,  p0(0, 2));
        EXPECT_EQ(0.0,  p0(0, 3));
        EXPECT_EQ(0.0,  p0(0, 4));

        EXPECT_EQ(0.6,  p0(1, 0));
        EXPECT_EQ(1.5,  p0(1, 1));
        EXPECT_EQ(0.0,  p0(1, 2));
        EXPECT_EQ(0.0,  p0(1, 3));
        EXPECT_EQ(0.0,  p0(1, 4));

        EXPECT_EQ(0.2,  p0(2, 0));
        EXPECT_EQ(0.9,  p0(2, 1));
        EXPECT_EQ(1.92, p0(2, 2));
        EXPECT_EQ(0.0,  p0(2, 3));
        EXPECT_EQ(0.0,  p0(2, 4));

        EXPECT_EQ(0.2,  p0(3, 0));
        EXPECT_EQ(0.3,  p0(3, 1));
        EXPECT_EQ(0.96, p0(3, 2));
        EXPECT_EQ(2.8,  p0(3, 3));
        EXPECT_EQ(1.06, p0(3, 4));

        EXPECT_EQ(0.2,  p0(4, 0));
        EXPECT_EQ(0.3,  p0(4, 1));
        EXPECT_EQ(0.32, p0(4, 2));
        EXPECT_EQ(1.2,  p0(4, 3));
        EXPECT_EQ(9.54, p0(4, 4));
    }
}

//---------------------------------------------------------------------------//

TEST_F(XS_Builder_Test, load_5GP1_P0_G13)
{
    builder.open_and_broadcast("xs5GP1.xml");
    EXPECT_EQ(5, builder.num_groups());
    EXPECT_EQ(1, builder.pn_order());

    const Vec_Str &mats = builder.materials();
    EXPECT_EQ(3, mats.size());
    EXPECT_EQ("mat 1", mats[0]);
    EXPECT_EQ("mat 4", mats[1]);
    EXPECT_EQ("mat 5", mats[2]);

    // read 2 materials
    Matid_Map map;
    map.insert(Matid_Map::value_type(1, std_string("mat 4")));
    map.insert(Matid_Map::value_type(2, std_string("mat 5")));
    map.complete();

    // build the cross sections
    builder.build(map, 0, 1, 3);
    RCP_XS xs = builder.get_xs();
    EXPECT_FALSE(xs.is_null());

    // test the xs
    EXPECT_EQ(0, xs->pn_order());
    EXPECT_EQ(3, xs->num_groups());
    EXPECT_EQ(2, xs->num_mat());

    EXPECT_TRUE(xs->has(1));
    EXPECT_TRUE(xs->has(2));

    // test bounds
    const auto &bounds = xs->bounds();
    EXPECT_EQ(4, bounds.length());
    EXPECT_EQ(1.0e6,  bounds(0));
    EXPECT_EQ(1.0e5,  bounds(1));
    EXPECT_EQ(1.0e4,  bounds(2));
    EXPECT_EQ(1.0e3,  bounds(3));

    // material 2
    {
        const Vector &sigt = xs->vector(2, XS_t::TOTAL);
        EXPECT_EQ(3, sigt.length());

        EXPECT_EQ(7.5, sigt(0));
        EXPECT_EQ(8, sigt(1));
        EXPECT_EQ(10, sigt(2));

        const Vector &sigf = xs->vector(2, XS_t::SIG_F);
        EXPECT_EQ(3, sigf.length());

        EXPECT_EQ(2.25, sigf(0));
        EXPECT_EQ(2.4, sigf(1));
        EXPECT_EQ(3.0, sigf(2));

        const Vector &nusigf = xs->vector(2, XS_t::NU_SIG_F);
        EXPECT_EQ(3, nusigf.length());

        EXPECT_EQ(5.4, nusigf(0));
        EXPECT_EQ(5.76, nusigf(1));
        EXPECT_EQ(7.2, nusigf(2));

        const Vector &chi = xs->vector(2, XS_t::CHI);
        EXPECT_EQ(3, chi.length());

        EXPECT_EQ(0.3, chi(0));
        EXPECT_EQ(0.1, chi(1));
        EXPECT_EQ(0.0, chi(2));

        const Matrix &p0 = xs->matrix(2, 0);
        EXPECT_EQ(3, p0.numRows());
        EXPECT_EQ(3, p0.numCols());

        EXPECT_EQ(1.5,  p0(0, 0));
        EXPECT_EQ(0.0,  p0(0, 1));
        EXPECT_EQ(0.0,  p0(0, 2));

        EXPECT_EQ(0.9,  p0(1, 0));
        EXPECT_EQ(1.92, p0(1, 1));
        EXPECT_EQ(0.0,  p0(1, 2));

        EXPECT_EQ(0.3,  p0(2, 0));
        EXPECT_EQ(0.96, p0(2, 1));
        EXPECT_EQ(2.8,  p0(2, 2));
    }
}

//---------------------------------------------------------------------------//
//                 end of tstXS_Builder.cc
//---------------------------------------------------------------------------//
