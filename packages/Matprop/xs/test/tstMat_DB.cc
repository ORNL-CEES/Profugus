//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Matprop/xs/test/tstMat_DB.cc
 * \author Thomas M. Evans
 * \date   Mon Feb 10 20:20:03 2014
 * \brief  Mat_DB unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../XS_Builder.hh"
#include "../Mat_DB.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class Mat_DB_Test : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture
    typedef profugus::XS_Builder  XS_Builder;
    typedef profugus::Mat_DB      Mat_DB;
    typedef Mat_DB::RCP_XS        RCP_XS;
    typedef XS_Builder::Matid_Map Matid_Map;
    typedef XS_Builder::Vec_Str   Vec_Str;
    typedef Mat_DB::Vec_Int       Vec_Int;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        // matid map
        Matid_Map map;

        // load cross sections
        XS_Builder builder;
        builder.open_and_broadcast("xs5GP1.xml");

        const Vec_Str &mats = builder.materials();
        int matid           = 1;
        for (Vec_Str::const_iterator m = mats.begin(); m != mats.end(); ++m)
        {
            map.insert(Matid_Map::value_type(matid, *m));
            ++matid;
        }
        map.complete();

        // build the cross sections
        builder.build(map);
        xs = builder.get_xs();

        // assign
        mat.set(xs, 4);
        mat.matid(0) = 1;
        mat.matid(1) = 3;
        mat.matid(2) = 2;
        mat.matid(3) = 1;
    }

  protected:
    // >>> Data that get re-initialized between tests

    RCP_XS xs;
    Mat_DB mat;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Mat_DB_Test, assign)
{
    mat.set(xs, 5);
    EXPECT_EQ(5, mat.num_cells());
    for (int i = 0; i < 5; ++i)
    {
        EXPECT_EQ(-1, mat.matid(i));
    }

    Vec_Int matids(3);
    matids[0] = 1;
    matids[1] = 2;
    matids[2] = 3;

    mat.assign(matids);
    EXPECT_EQ(3, mat.num_cells());
    EXPECT_EQ(1, mat.matid(0));
    EXPECT_EQ(2, mat.matid(1));
    EXPECT_EQ(3, mat.matid(2));

    mat.set(xs);
    EXPECT_EQ(0, mat.num_cells());
}

//---------------------------------------------------------------------------//

TEST_F(Mat_DB_Test, access)
{
    EXPECT_EQ(1, mat.matid(0));
    EXPECT_EQ(3, mat.matid(1));
    EXPECT_EQ(2, mat.matid(2));
    EXPECT_EQ(1, mat.matid(3));

    const Vec_Int &matids = mat.matids();
    EXPECT_EQ(4, matids.size());
    EXPECT_EQ(1, matids[0]);
    EXPECT_EQ(3, matids[1]);
    EXPECT_EQ(2, matids[2]);
    EXPECT_EQ(1, matids[3]);

    EXPECT_EQ(4, mat.num_cells());
}

//---------------------------------------------------------------------------//
//                 end of tstMat_DB.cc
//---------------------------------------------------------------------------//
