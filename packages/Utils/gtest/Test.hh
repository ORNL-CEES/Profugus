//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/gtest/Test.hh
 * \author Seth R Johnson
 * \date   Fri Mar 11 07:43:01 2016
 * \brief  Test class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_gtest_Test_hh
#define Utils_gtest_Test_hh

#include <string>
#include "gtest.h"

namespace profugus
{

//===========================================================================//
/*!
 * \class Test
 * \brief MPI-aware google test case.
 *
 * This provides access to "node" and "nodes" data attributes, and it provides
 * a convenience function for making a test-unique filename.
 *
 * \code
    TEST_F(Things, stuff)
    {
        std::string my_filename = get_filename(".h5");
        // In serial mode, returns "things-stuff.h5"
        // In parallel mode, returns e.g. "things-stuff-np2.h5"
    }
 * \endcode
 */
/*!
 * \example gtest/test/tstTestCase.cc
 *
 * Test of Test.
 */
//===========================================================================//

class Test : public ::testing::Test
{
  protected:
    // >>> PROTECED DATA

    int node;
    int nodes;

  private:

    // >>> PRIVATE DATA
    int d_filename_counter;

  public:

    // Constructor
    Test();

    // Generate test-unique filename
    std::string make_unique_filename(const char* ext = "");

};

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
#endif // Utils_gtest_Test_hh

//---------------------------------------------------------------------------//
// end of Utils/gtest/Test.hh
//---------------------------------------------------------------------------//
