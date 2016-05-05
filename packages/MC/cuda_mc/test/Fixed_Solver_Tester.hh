//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fixed_Solver_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Fixed_Solver_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Fixed_Solver_Tester_hh
#define cuda_mc_Fixed_Solver_Tester_hh

//===========================================================================//
/*!
 * \class Fixed_Solver_Tester
 * \brief Wrapper to kernel launches to test Fixed_Solver class.
 */
//===========================================================================//

class Fixed_Solver_Tester
{
  public:

      static void test_transport(int num_groups);
};

#endif // cuda_mc_Fixed_Solver_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Fixed_Solver_Tester.hh
//---------------------------------------------------------------------------//
