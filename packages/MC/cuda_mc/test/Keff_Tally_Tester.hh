//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Keff_Tally_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Keff_Tally_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Keff_Tally_Tester_hh
#define cuda_mc_Keff_Tally_Tester_hh

//===========================================================================//
/*!
 * \class Keff_Tally_Tester
 * \brief Wrapper to kernel launches to test Keff_Tally class.
 */
//===========================================================================//

class Keff_Tally_Tester
{
  public:

      static void test_tally(int num_groups);
};

#endif // cuda_mc_Keff_Tally_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Keff_Tally_Tester.hh
//---------------------------------------------------------------------------//
