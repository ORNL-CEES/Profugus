//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/KCode_Solver_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  KCode_Solver_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_KCode_Solver_Tester_hh
#define cuda_mc_KCode_Solver_Tester_hh

//===========================================================================//
/*!
 * \class KCode_Solver_Tester
 * \brief Wrapper to kernel launches to test KCode_Solver class.
 */
//===========================================================================//

class KCode_Solver_Tester
{
  public:

      static void test_mesh(int num_groups);
      static void test_rtk();
};

#endif // cuda_mc_KCode_Solver_Tester_hh

//---------------------------------------------------------------------------//
//                 end of KCode_Solver_Tester.hh
//---------------------------------------------------------------------------//
