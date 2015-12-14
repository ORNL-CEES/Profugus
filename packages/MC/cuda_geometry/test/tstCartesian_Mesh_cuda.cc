//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/test/tstCartesian_Mesh.cc
 * \author hu4
 * \date   Mon Dec 14 13:16:35 2015
 * \brief  Cartesian_Mesh class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Cartesian_Mesh_Tester.hh"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class CartesianMeshTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef std::vector<int>                     Vec_Int;
    typedef std::vector<double>                  Vec_Dbl;
    typedef cuda_profugus::Cartesian_Mesh_Tester Tester;

  protected:
    void SetUp()
    {
        Vec_Dbl x_edges = {0.0, 0.1, 0.6, 0.9, 1.0};
        Vec_Dbl y_edges = {-1.0, -0.6, 0.0};
        Vec_Dbl z_edges = {2.0, 2.6, 3.4, 4.0};

        d_tester = std::make_shared<Tester>(x_edges,y_edges,z_edges);
    }

  protected:
    // >>> DATA
    std::shared_ptr<Tester> d_tester;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(CartesianMeshTest, Volumes)
{

    Vec_Int cells = {4, 1, 22, 11};
    Vec_Dbl volumes(cells.size());

    d_tester->compute_volumes(cells,volumes);

    std::cout << "Computed volumes: ";
    for( auto x : volumes )
        std::cout << x << " ";
    std::cout << std::endl;

    Vec_Dbl expected_volumes = {0.1 * 0.6 * 0.6,
                                0.5 * 0.4 * 0.6,
                                0.3 * 0.6 * 0.6,
                                0.1 * 0.4 * 0.8};

    EXPECT_VEC_SOFTEQ( expected_volumes, volumes, 1e-12 );
}

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/tstCartesian_Mesh.cc
//---------------------------------------------------------------------------//
