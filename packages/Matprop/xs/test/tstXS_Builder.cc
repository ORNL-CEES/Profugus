//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/test/tstXS_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Feb 05 20:07:26 2014
 * \brief  XS_Builder unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../XS_Builder.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class XS_Builder_Test : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture
    typedef profugus::XS_Builder XS_Builder;
    typedef XS_Builder::XS_t     XS_t;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();
    }

  protected:
    // >>> Data that get re-initialized between tests

    int node, nodes;
    XS_Builder builder;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(XS_Builder_Test, load_p0)
{
    builder.open_and_broadcast("xs3GP0.xml");
}

//---------------------------------------------------------------------------//
//                 end of tstXS_Builder.cc
//---------------------------------------------------------------------------//
