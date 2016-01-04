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

#include <vector>
#include <memory>

#include "cuda_utils/Definitions.hh"
#include "geometry/Definitions.hh"

namespace cuda_profugus
{

class Mesh_Geometry;

//===========================================================================//
/*!
 * \class Mesh_Geometry_Tester
 * \brief Helper class for testing Mesh_Geometry kernels.
 *
 */
//===========================================================================//

class Mesh_Geometry_Tester
{
  public:

      typedef profugus::geometry::cell_type  cell_type;
      typedef profugus::geometry::matid_type matid_type;
      typedef std::vector<double>            Vec_Dbl;
      typedef std::vector<int>               Vec_Int;
      typedef std::vector<cell_type>         Vec_Cell_Type;
      typedef std::vector<matid_type>        Vec_Matid_Type;
      typedef cuda::Space_Vector             Point;
      typedef std::vector<Point>             Vec_Point;

      Mesh_Geometry_Tester( const Vec_Dbl &x_edges,
                             const Vec_Dbl &y_edges,
                             const Vec_Dbl &z_edges );

      void compute_volumes(const Vec_Cell_Type &cells,
                                 Vec_Dbl       &volumes) const;

      void compute_matids(const Vec_Matid_Type &all_matids,
                          const Vec_Point      &points,
                                Vec_Matid_Type &matids) const;

  private:

      std::shared_ptr<Mesh_Geometry> d_mesh;

};

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

#endif // MC_cuda_geometry_test_Mesh_Geometry_Tester_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/Mesh_Geometry_Tester.hh
//---------------------------------------------------------------------------//
