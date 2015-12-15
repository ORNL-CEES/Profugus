//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/test/tstCartesian_Mesh.cc
 * \author Steven Hamilton
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

TEST_F(CartesianMeshTest, Index)
{

    Vec_Int ii = {0, 1, 2, 3};
    Vec_Int jj = {1, 0, 1, 0};
    Vec_Int kk = {0, 0, 2, 1};
    Vec_Int cells(ii.size());

    d_tester->compute_indices(ii,jj,kk,cells);

    std::cout << "Computed indices: ";
    for( auto x : cells )
        std::cout << x << " ";
    std::cout << std::endl;
    Vec_Int expected_cells = {4, 1, 22, 11};

    EXPECT_VEC_EQ( expected_cells, cells );
}

TEST_F(CartesianMeshTest, Cardinal)
{

    Vec_Int cells = {4, 1, 22, 11};
    Vec_Int ii(cells.size());
    Vec_Int jj(cells.size());
    Vec_Int kk(cells.size());

    d_tester->compute_cardinals(cells,ii,jj,kk);

    std::cout << "Computed cardinals: " << std::endl;
    for( int cell = 0; cell < cells.size(); ++cell )
    {
        std::cout << cells[cell] << ": " << ii[cell] << " " << jj[cell]
            << " " << kk[cell] << std::endl;
    }

    Vec_Int expected_ii = {0, 1, 2, 3};
    Vec_Int expected_jj = {1, 0, 1, 0};
    Vec_Int expected_kk = {0, 0, 2, 1};

    EXPECT_VEC_EQ( expected_ii, ii );
    EXPECT_VEC_EQ( expected_jj, jj );
    EXPECT_VEC_EQ( expected_kk, kk );
}

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

    EXPECT_VEC_SOFT_EQ( expected_volumes, volumes );
}

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/tstCartesian_Mesh.cc
//---------------------------------------------------------------------------//
