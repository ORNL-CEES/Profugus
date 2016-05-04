//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Cell_Tally_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Cell_Tally_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Cell_Tally_Tester_hh
#define cuda_mc_Cell_Tally_Tester_hh

//===========================================================================//
/*!
 * \class Cell_Tally_Tester
 * \brief Wrapper to kernel launches to test Cell_Tally class.
 */
//===========================================================================//

class Cell_Tally_Tester
{
  public:


      static void test_tally();
};

#endif // cuda_mc_Cell_Tally_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Cell_Tally_Tester.hh
//---------------------------------------------------------------------------//
