//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/gtest/Test.cc
 * \author Seth R Johnson
 * \date   Fri Mar 11 07:43:01 2016
 * \brief  Test class definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Test.hh"

#include <sstream>
#include "Utils/harness/DBC.hh"
#include "Utils/comm/global.hh"
#include "Utils/utils/String_Functions.hh"

namespace profugus
{
//---------------------------------------------------------------------------//
/*!
 * \brief Set node and nodes
 */
Test::Test()
    : node(profugus::node())
    , nodes(profugus::nodes())
    , d_filename_counter(0)
{
    /* * */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Generate test-unique filename.
 */
std::string Test::make_unique_filename(const char* ext)
{
    REQUIRE(ext);

    // Get filename based on unit test name
    const ::testing::TestInfo* const test_info =
        ::testing::UnitTest::GetInstance()->current_test_info();
    CHECK(test_info);

    // Convert test case to lowercase
    std::string case_name = profugus::lower(test_info->test_case_name());

    // Delete "Test" if present
    auto pos = case_name.find("test");
    if (pos != std::string::npos)
    {
        case_name.replace(pos, 4, "", 0);
    }

    // Delete leading underscore
    if (!case_name.empty() && case_name.front() == '_')
    {
        case_name.erase(case_name.begin(), case_name.begin() + 1);
    }

    std::ostringstream os;
    os << case_name << '-' << test_info->name();

    if (d_filename_counter)
    {
        os << '-' << d_filename_counter;
    }
    ++d_filename_counter;

#ifdef USE_MPI
    os << "-np" << nodes;
#endif
    os << ext;

    return os.str();
}

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
// end of Utils/gtest/Test.cc
//---------------------------------------------------------------------------//
