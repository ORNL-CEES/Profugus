//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Box_Shape_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 14:57:16 2016
 * \brief  Box_Shape_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Box_Shape_Tester_hh
#define cuda_mc_Box_Shape_Tester_hh

#include <vector>
#include "cuda_utils/Definitions.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Box_Shape_Tester
 * \brief Wrapper to kernel launch to test Box_Shape class.
 */
//===========================================================================//

class Box_Shape_Tester
{

  public:

      typedef std::vector<int>                Vec_Int;
      typedef std::vector<double>             Vec_Dbl;
      typedef std::vector<cuda::Space_Vector> Vec_Space_Vec;

      static void test_inside( const Vec_Dbl       &box_bounds,
                               const Vec_Space_Vec &pts,
                                     Vec_Int       &inside );
      static void test_sample( const Vec_Dbl       &box_bounds,
                                     Vec_Space_Vec &pts );
};

} // end namespace cuda_mc

#endif // cuda_mc_Box_Shape_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Box_Shape_Tester.hh
//---------------------------------------------------------------------------//
