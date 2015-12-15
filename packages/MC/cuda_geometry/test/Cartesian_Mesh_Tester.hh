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

#include <vector>
#include <memory>

#include "geometry/Definitions.hh"

namespace cuda_profugus
{

class Cartesian_Mesh;

//===========================================================================//
/*!
 * \class Cartesian_Mesh_Tester
 * \brief Helper class for testing Cartesian_Mesh kernels.
 *
 */
//===========================================================================//

class Cartesian_Mesh_Tester
{
  public:

      typedef profugus::geometry::cell_type cell_type;
      typedef std::vector<double>           Vec_Dbl;
      typedef std::vector<int>              Vec_Int;
      typedef std::vector<cell_type>        Vec_Cell_Type;

      Cartesian_Mesh_Tester( const Vec_Dbl &x_edges,
                             const Vec_Dbl &y_edges,
                             const Vec_Dbl &z_edges );

      void compute_indices(const Vec_Int       &ii,
                           const Vec_Int       &jj,
                           const Vec_Int       &kk,
                                 Vec_Cell_Type &cells) const;

      void compute_cardinals(const Vec_Cell_Type &cells,
                                   Vec_Int       &ii,
                                   Vec_Int       &jj,
                                   Vec_Int       &kk ) const;

      void compute_volumes(const Vec_Cell_Type &cells,
                                 Vec_Dbl       &volumes) const;

  private:

      std::shared_ptr<Cartesian_Mesh> d_mesh;

};

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

#endif // MC_cuda_geometry_test_Cartesian_Mesh_Tester_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/Cartesian_Mesh_Tester.hh
//---------------------------------------------------------------------------//
