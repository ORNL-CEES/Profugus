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
#include "geometry/Definitions.hh"

TEST(CartesianMesh, index)
{
    Cartesian_Mesh_Tester::test_index();
}

TEST(CartesianMesh, cardinal)
{
    Cartesian_Mesh_Tester::test_cardinal();
}

TEST(CartesianMesh, volumes)
{
    Cartesian_Mesh_Tester::test_volume();
}

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/tstCartesian_Mesh.cc
//---------------------------------------------------------------------------//
