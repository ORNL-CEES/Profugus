//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Physics_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Physics_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Physics_Tester_hh
#define cuda_mc_Physics_Tester_hh


//===========================================================================//
/*!
 * \class Physics_Tester
 * \brief Wrapper to kernel launches to test Physics class.
 */
//===========================================================================//

class Physics_Tester
{
  public:

      static void test_total();
      static void test_collide();
};

#endif // cuda_mc_Physics_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Physics_Tester.hh
//---------------------------------------------------------------------------//
