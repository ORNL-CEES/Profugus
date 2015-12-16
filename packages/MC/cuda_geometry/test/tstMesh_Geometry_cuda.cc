//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/test/tstMesh_Geometry.cc
 * \author Steven Hamilton
 * \date   Mon Dec 14 13:16:35 2015
 * \brief  Mesh_Geometry class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Mesh_Geometry_Tester.hh"

#include "gtest/utils_gtest.hh"
#include "geometry/Definitions.hh"
#include "cuda_utils/Definitions.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class MeshGeometryTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::geometry::cell_type        cell_type;
    typedef profugus::geometry::matid_type       matid_type;
    typedef std::vector<int>                     Vec_Int;
    typedef std::vector<cell_type>               Vec_Cell_Type;
    typedef std::vector<matid_type>              Vec_Matid_Type;
    typedef std::vector<double>                  Vec_Dbl;
    typedef cuda::Space_Vector                   Point;
    typedef std::vector<Point>                   Vec_Point;
    typedef cuda_profugus::Mesh_Geometry_Tester  Tester;

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

TEST_F(MeshGeometryTest, Volumes)
{

    Vec_Cell_Type cells = {4, 1, 22, 11};
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

TEST_F(MeshGeometryTest, Matids)
{
    Vec_Matid_Type all_matids = {1, 3, 2, 0,
                                 3, 1, 4, 1,
                                 2, 5, 2, 1,
                                 0, 1, 2, 3,
                                 1, 2, 3, 4,
                                 2, 3, 4, 5};

    Vec_Point points = {{0.7,  -0.9,  2.1},
                        {0.5,  -0.5,  2.5},
                        {0.99, -0.01, 3.99},
                        {0.05, -0.8,  2.4}};
    Vec_Matid_Type cell_matids(points.size());

    d_tester->compute_matids(all_matids,points,cell_matids);

    std::cout << "Computed matids: ";
    for( auto x : cell_matids )
        std::cout << x << " ";
    std::cout << std::endl;

    Vec_Matid_Type expected_matids = {2, 1, 5, 1};

    EXPECT_VEC_EQ( expected_matids, cell_matids);
}


//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/tstMesh_Geometry.cc
//---------------------------------------------------------------------------//
