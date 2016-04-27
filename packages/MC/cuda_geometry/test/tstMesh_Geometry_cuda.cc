//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/test/tstMesh_Geometry.cc
 * \author Steven Hamilton
 * \date   Mon Dec 14 13:16:35 2015
 * \brief  Mesh_Geometry class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include "Mesh_Geometry_Tester.hh"

TEST(MeshGeometry, volume)
{
    Mesh_Geometry_Tester::test_volume();
}

TEST(MeshGeometry, matid)
{
    Mesh_Geometry_Tester::test_matid();
}

TEST(MeshGeometry, dist_to_bdry)
{
    Mesh_Geometry_Tester::test_dist_to_bdry();
}

TEST(MeshGeometry, move_to_surf)
{
    Mesh_Geometry_Tester::test_move_to_surf();
}

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/tstMesh_Geometry.cc
//---------------------------------------------------------------------------//
