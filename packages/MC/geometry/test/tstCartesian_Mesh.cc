//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/geometry/test/tstCartesian_Mesh.cc
 * \author Thomas M. Evans
 * \date   Mon Jul 21 17:24:52 2014
 * \brief  Cartesian_Mesh test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Cartesian_Mesh.hh"

#include "gtest/utils_gtest.hh"

using profugus::Cartesian_Mesh;
using def::I; using def::J; using def::K;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class CartesianMeshTest : public ::testing::Test
{
  protected:
    typedef Cartesian_Mesh       Mesh_t;
    typedef Mesh_t::Vec_Dbl      Vec_Dbl;
    typedef Mesh_t::size_type    size_type;
    typedef Mesh_t::Space_Vector Space_Vector;
    typedef Mesh_t::Dim_Vector   Dim_Vector;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        x.resize(5);
        x[0] = 0.0; x[1] = 0.10; x[2] = 0.25; x[3] = 0.30; x[4] = 0.42;

        y.resize(4);
        y[0] = 0.0; y[1] = 0.20; y[2] = 0.40; y[3] = 0.50;

        z.resize(4);
        z[0] = -0.1; z[1] = 0.0; z[2] = 0.15; z[3] = 0.50;
    }

  protected:
    // >>> Data that get re-initialized between tests
    Vec_Dbl x;
    Vec_Dbl y;
    Vec_Dbl z;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(CartesianMeshTest, accessors)
{
    Mesh_t mesh(x, y, z);

    EXPECT_EQ(36, mesh.num_cells());
    EXPECT_EQ(4 , mesh.num_cells_along(I));
    EXPECT_EQ(3 , mesh.num_cells_along(J));
    EXPECT_EQ(3 , mesh.num_cells_along(K));

    EXPECT_EQ(3 , mesh.dimension());

    EXPECT_SOFT_EQ( 0. , mesh.low_corner(I));
    EXPECT_SOFT_EQ( 0. , mesh.low_corner(J));
    EXPECT_SOFT_EQ(-0.1, mesh.low_corner(K));

    EXPECT_SOFT_EQ(0.42, mesh.high_corner(I));
    EXPECT_SOFT_EQ(0.50, mesh.high_corner(J));
    EXPECT_SOFT_EQ(0.50, mesh.high_corner(K));

    EXPECT_EQ(Dim_Vector(4,3,3), mesh.extents());

    // Test index/cardinal
    for (size_t i = 0; i < mesh.num_cells(); ++i)
    {
        Dim_Vector ijk = mesh.cardinal(i);

        size_type  cell;
        bool success = mesh.index(ijk[I], ijk[J], ijk[K], cell);
        ASSERT_TRUE(success);
        EXPECT_EQ(i, cell);
        // Alternate index call
        EXPECT_EQ(i, mesh.index(ijk[I], ijk[J], ijk[K]));
    }

    // Test volume
    EXPECT_SOFTEQ(0.002 , mesh.volume(0, 0, 0) , 1.e-12);

    double expected_vol[] = {0.002 , 0.003 , 0.001 , 0.0024, 0.002 , 0.003 ,
        0.001 , 0.0024, 0.001 , 0.0015, 0.0005, 0.0012, 0.003 , 0.0045, 0.0015,
        0.0036};

    for (size_t i = 0; i < 16; ++i)
    {
        EXPECT_SOFTEQ(expected_vol[i], mesh.volume(i), 1.e-12)
            << "Volume mismatch in cell " << i;
    }
}

//---------------------------------------------------------------------------//

TEST_F(CartesianMeshTest, find)
{
    Mesh_t mesh(x, y, z);

    // Check outside position
    Space_Vector r(-1., 1., 0.1);

    Dim_Vector ijk;
    mesh.find_upper(r, ijk);
    EXPECT_EQ(-1 , ijk[I]);
    EXPECT_EQ( 3 , ijk[J]);
    EXPECT_EQ( 1 , ijk[K]);

    bool success;
    success = mesh.find(r, ijk);
    EXPECT_FALSE(success);

    // Check inside position
    r = Space_Vector(0.2, 0.45, 0.4);
    mesh.find_upper(r, ijk);
    EXPECT_EQ(1 , ijk[I]);
    EXPECT_EQ(2 , ijk[J]);
    EXPECT_EQ(2 , ijk[K]);

    ijk = Dim_Vector(-1,-1,-1);
    success = mesh.find(r, ijk);
    EXPECT_TRUE(success);
    EXPECT_EQ(1 , ijk[I]);
    EXPECT_EQ(2 , ijk[J]);
    EXPECT_EQ(2 , ijk[K]);
}

//---------------------------------------------------------------------------//
//                 end of tstCartesian_Mesh.cc
//---------------------------------------------------------------------------//
