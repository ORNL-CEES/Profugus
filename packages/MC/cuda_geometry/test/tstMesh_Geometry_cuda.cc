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

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/tstMesh_Geometry.cc
//---------------------------------------------------------------------------//
