//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Source_Transporter_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Source_Transporter_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Source_Transporter_Tester_hh
#define cuda_mc_Source_Transporter_Tester_hh

//===========================================================================//
/*!
 * \class Source_Transporter_Tester
 * \brief Wrapper to kernel launches to test Source_Transporter class.
 */
//===========================================================================//

class Source_Transporter_Tester
{
  public:

      static void test_transport(int num_groups);
};

#endif // cuda_mc_Source_Transporter_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Source_Transporter_Tester.hh
//---------------------------------------------------------------------------//
