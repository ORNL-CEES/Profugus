//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Uniform_Source_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 14:57:16 2016
 * \brief  Uniform_Source_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Uniform_Source_Tester_hh
#define cuda_mc_Uniform_Source_Tester_hh

#include <vector>
#include "cuda_utils/Definitions.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Uniform_Source_Tester
 * \brief Wrapper to kernel launch to test Uniform_Source class.
 */
//===========================================================================//

class Uniform_Source_Tester
{

  public:

      typedef std::vector<double>             Vec_Dbl;
      typedef std::vector<cuda::Space_Vector> Vec_Space_Vec;

      static void test_source( const Vec_Dbl    &geom_bounds,
                               const Vec_Dbl    &src_bounds,
                                     Vec_Dbl    &x_loc,
                                     Vec_Dbl    &y_loc,
                                     Vec_Dbl    &z_loc,
                                     Vec_Dbl    &mu,
                                     Vec_Dbl    &eta,
                                     Vec_Dbl    &xi);

      static void test_host_api( const Vec_Dbl    &geom_bounds,
                                 const Vec_Dbl    &src_bounds,
                                       Vec_Dbl    &x_loc,
                                       Vec_Dbl    &y_loc,
                                       Vec_Dbl    &z_loc,
                                       Vec_Dbl    &mu,
                                       Vec_Dbl    &eta,
                                       Vec_Dbl    &xi);
};

} // end namespace cuda_mc

#endif // cuda_mc_Uniform_Source_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Uniform_Source_Tester.hh
//---------------------------------------------------------------------------//
