//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/test/Cartesian_Mesh_Tester.hh
 * \author Steven Hamilton
 * \date   Mon Dec 14 13:28:26 2015
 * \brief  Cartesian_Mesh_Tester class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_geometry_test_Cartesian_Mesh_Tester_hh
#define MC_cuda_geometry_test_Cartesian_Mesh_Tester_hh

//===========================================================================//
/*!
 * \class Cartesian_Mesh_Tester
 */
//===========================================================================//

class Cartesian_Mesh_Tester
{
  public:

      static void test_index();
      static void test_cardinal();
      static void test_volume();
      static void test_find_upper();
};

#endif // MC_cuda_geometry_test_Cartesian_Mesh_Tester_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/Cartesian_Mesh_Tester.hh
//---------------------------------------------------------------------------//
