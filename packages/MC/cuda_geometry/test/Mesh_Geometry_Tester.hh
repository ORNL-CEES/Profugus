//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/test/Mesh_Geometry_Tester.hh
 * \author Steven Hamilton
 * \date   Mon Dec 14 13:28:26 2015
 * \brief  Mesh_Geometry_Tester class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_geometry_test_Mesh_Geometry_Tester_hh
#define MC_cuda_geometry_test_Mesh_Geometry_Tester_hh

//===========================================================================//
/*!
 * \class Mesh_Geometry_Tester
 */
//===========================================================================//

class Mesh_Geometry_Tester
{
  public:

      static void test_volume();
      static void test_matid();
      //static void test_dist_to_bdry();
};

#endif // MC_cuda_geometry_test_Mesh_Geometry_Tester_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/Mesh_Geometry_Tester.hh
//---------------------------------------------------------------------------//
