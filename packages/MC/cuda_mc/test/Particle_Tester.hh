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

#include <vector>

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Particle_Tester
 * \brief Wrapper to kernel launch to test Particle class.
 */
//===========================================================================//

class Particle_Tester
{

  public:

      typedef std::vector<double>   Vec_Dbl;
      typedef std::vector<int>      Vec_Int;

      static void test_randoms( Vec_Dbl &randoms );
      static void test_groups( const Vec_Int &groups_in, Vec_Int &groups_out);
      static void test_matids( const Vec_Int &matids_in, Vec_Int &matids_out);

};

} // end namespace cuda_mc

#endif // cuda_mc_Particle_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Particle_Tester.hh
//---------------------------------------------------------------------------//
