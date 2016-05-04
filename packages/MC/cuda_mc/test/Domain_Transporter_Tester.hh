//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Domain_Transporter_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Domain_Transporter_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Domain_Transporter_Tester_hh
#define cuda_mc_Domain_Transporter_Tester_hh

//===========================================================================//
/*!
 * \class Domain_Transporter_Tester
 * \brief Wrapper to kernel launches to test Domain_Transporter class.
 */
//===========================================================================//

class Domain_Transporter_Tester
{
  public:

      static void test_transport(int num_groups);
};

#endif // cuda_mc_Domain_Transporter_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter_Tester.hh
//---------------------------------------------------------------------------//
