//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Particle_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 14:57:16 2016
 * \brief  Particle_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Particle_Tester_hh
#define cuda_mc_Particle_Tester_hh

//===========================================================================//
/*!
 * \class Particle_Tester
 * \brief Wrapper to kernel launch to test Particle class.
 */
//===========================================================================//

class Particle_Tester
{

  public:

      static void test_randoms();
      static void test_groups();
      static void test_matids();

};

#endif // cuda_mc_Particle_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Particle_Tester.hh
//---------------------------------------------------------------------------//
