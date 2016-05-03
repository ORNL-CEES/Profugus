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

//===========================================================================//
/*!
 * \class Box_Shape_Tester
 * \brief Wrapper to kernel launch to test Box_Shape class.
 */
//===========================================================================//

class Box_Shape_Tester
{

  public:

      static void test_inside();
      static void test_sample();
};

#endif // cuda_mc_Box_Shape_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Box_Shape_Tester.hh
//---------------------------------------------------------------------------//
