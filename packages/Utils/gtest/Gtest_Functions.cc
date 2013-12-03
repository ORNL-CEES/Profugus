//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   gtest/Gtest_Functions.cc
 * \author Seth R Johnson
 * \date   Tue Apr 02 11:41:41 2013
 * \brief  Gtest_Functions member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Gtest_Functions.hh"

#include <gtest/config.h>
#include <vector>
#include <iomanip>
#include <sstream>
#include <limits>
#include <cmath>
#include <cstddef>

#include "harness/Soft_Equivalence.hh"
#include "harness/DBC.hh"
#include "harness/Warnings.hh"
#include "release/Release.hh"
#include "comm/global.hh"

#include "gtest.h"
#include "gtest-internals.hh"

//---------------------------------------------------------------------------//
// ANONYMOUS NAMESPACE HELPER FUNCTIONS
//---------------------------------------------------------------------------//
namespace
{
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the number of base-10 digits a nonnegative number has
 */
unsigned int calc_num_digits(std::size_t number)
{
    unsigned int num_digits = 0;
    do {
        number /= 10;
        ++num_digits;
    } while (number > 0);

    Ensure(num_digits > 0);
    return num_digits;
}
}
//---------------------------------------------------------------------------//

namespace nemesis
{
//---------------------------------------------------------------------------//
// INTERNAL-USE CLAASES
//---------------------------------------------------------------------------//
class ParallelHandler : public ::testing::EmptyTestEventListener
{
  public:
    // This will only be called after MPI_Init, so we can access these data
    ParallelHandler()
      : d_node(nemesis::node())
      , d_num_nodes(nemesis::nodes())
    {
        /* * */
    }

  public:
    virtual void OnTestStart(const ::testing::TestInfo& test_info)
    {
        // Barrier before starting the run
        nemesis::global_barrier();
    }

    virtual void OnTestEnd(const ::testing::TestInfo& test_info);

  private:
    int d_node;
    int d_num_nodes;
};

//---------------------------------------------------------------------------//
void ParallelHandler::OnTestEnd(const ::testing::TestInfo& test_info)
{
    using ::testing::internal::ColoredPrintf;
    using ::testing::internal::COLOR_YELLOW;

    // Barrier after finishing the test part
    nemesis::global_barrier();

    // Print warnings
    if (!NEMESIS_WARNINGS.empty())
    {
        ColoredPrintf(COLOR_YELLOW,
                "%d WARNING%s noted on node %d in subtest '%s'\n",
                  NEMESIS_WARNINGS.num_warnings(),
                  NEMESIS_WARNINGS.num_warnings() > 1 ? "S" : "",
                  d_node,
                  test_info.name());
    }

    if (d_node == 0)
    {
        // Print warnings
        while (!NEMESIS_WARNINGS.empty())
        {
            ColoredPrintf(COLOR_YELLOW, "*** ");
            std::cout << NEMESIS_WARNINGS.pop() << std::endl;
        }
    }
    else
    {
        // Suppress warnings on other processors
        NEMESIS_WARNINGS.clear();
    }

    nemesis::global_barrier();
}

//---------------------------------------------------------------------------//
// NEMESIS FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \fn gtest_main
 * \brief Implementation for the harness' "main" function.
 *
 * This should be called by the main() function in each Google-based unit test.
 * It handles MPI initialization, the running and "finalizing" of the unit test
 * data (to ensure that if only one processor on a multi-processor run fails,
 * the overall result is a failure), etc. It also prints warnings at the end of
 * each test.
 */
int gtest_main(int argc, char *argv[])
{
    // Initialize MPI
    nemesis::initialize(argc, argv);
    const int node      = nemesis::node();
    const int num_nodes = nemesis::nodes();

    if (node == 0)
    {
        std::cout << "Using " << num_nodes << " processors" << std::endl
            << "Exnihilo " << nemesis::release::long_version() << std::endl;
    }
    else
    {
        // Disable color on non-main processor for visibility
        ::testing::GTEST_FLAG(color) = "no";
    }

    // Initialize google test
    ::testing::InitGoogleTest(&argc, argv);

    // Gets hold of the event listener list.
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();
    // Adds a listener to the end.  Google Test takes the ownership.
    listeners.Append(new ParallelHandler());

    // Run everything
    int failed = RUN_ALL_TESTS();

    // Accumulate the result so that all processors will have the same result
    nemesis::global_sum(failed);

    // Finish MPI
    nemesis::global_barrier();
    nemesis::finalize();

    // Print final results
    if (node == 0)
    {
        // Print warnings
        while (!NEMESIS_WARNINGS.empty())
        {
            std::cout << NEMESIS_WARNINGS.pop() << std::endl;
        }

        if (argc)
            std::cout << "In " << argv[0] << ", ";
        std::cout << "overall test result: "
            << (failed ? "FAILED" : "PASSED")
            << std::endl;
    }

    // Return 1 if any failure, 0 if all success
    return (failed ? 1 : 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Print the "skip" message from the skip macro
 *
 * This uses an internal function exposed by our small "gtest-internals.hh"
 * header. Note that if we update our version of google test, we'll probably
 * have to change this implementation a bit.
 */
void print_skip_message()
{
    ::testing::internal::ColoredPrintf(
            ::testing::internal::COLOR_YELLOW,
            "[   SKIP   ] ");
}

//---------------------------------------------------------------------------//
/*
 * \brief Custom error mesages for relative error soft equiavelence
 */
::testing::AssertionResult IsSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        const char* eps_expr,
        double expected,
        double actual,
        double eps)
{
    if (nemesis::soft_equiv(actual, expected, eps))
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        ::testing::AssertionResult failure = ::testing::AssertionFailure();

        failure << "Value of: " << actual_expr << "\n"
                << "  Actual: " << actual << "\n"
                << "Expected: " << expected_expr << "\n"
                << "Which is: " << expected << "\n";

        if (std::fabs(expected) < 1.e-14)
        {
            // Avoid divide by zero errors
            failure << "(Absolute error " << actual - expected;
        }
        else
        {
            failure << "(Relative error " << (actual - expected) / expected;
        }
        failure << " exceeds tolerance " << eps_expr << ")";

        return failure;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Custom vector comparison with soft equiavelence
 */
::testing::AssertionResult IsIterSoftEquiv(
        const char* expected_expr,
        const char* actual_expr,
        const char* eps_expr,
        const double* expected, const double* expected_end,
        const double* actual,   const double* actual_end,
        double eps)
{
    typedef std::size_t size_type;

    size_type expected_size = static_cast<size_type>(expected_end - expected);
    size_type actual_size   = static_cast<size_type>(actual_end - actual);

    // First, check that the sizes are equal
    if (expected_size != actual_size)
    {
        ::testing::AssertionResult failure = ::testing::AssertionFailure();

        failure << " Size of: " << actual_expr << "\n"
                << "  Actual: " << actual_size << "\n"
                << "Expected: " << expected_expr << ".size()\n"
                << "Which is: " << expected_size << "\n";
        return failure;
    }

    typedef std::vector<size_type> Vec_Size_Type;

    // Keep track of what elements failed
    Vec_Size_Type failed_indices;

    // Now loop through elements of the vectors and check for soft equivalence
    for (size_type i = 0; i < expected_size; ++i)
    {
        if (!nemesis::soft_equiv(actual[i], expected[i], eps))
        {
            failed_indices.push_back(i);
        }
    }

    if (failed_indices.empty())
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        using std::setw;
        using std::setprecision;

        ::testing::AssertionResult failure = ::testing::AssertionFailure();

        failure << "Values in: " << actual_expr << "\n"
                << " Expected: " << expected_expr << "\n"
                << failed_indices.size() << " of "
                << expected_size << " elements differ more than "
                << "given tolerance " << eps_expr << ":\n";

        // Only print the first 40 failures
        Vec_Size_Type::const_iterator end_indices = failed_indices.end();
        if (failed_indices.size() > 40)
        {
            failure << "(Truncating to first 40 failed values)\n";
            end_indices = failed_indices.begin() + 40;
        }

        // Calculate how many digits we need to space out
        unsigned int num_digits = 0;
        {
            size_type temp = failed_indices.back();
            do {
                temp /= 10;
                ++num_digits;
            } while (temp > 0);
        }

        // Construct our own stringstream because google test ignores setw
        std::ostringstream failure_stream;
        failure_stream << setprecision(std::numeric_limits<double>::digits10);

        double error = -1.;

        // Try to use user-given expressions for headers, but fall back if the
        // column length is exceeded
        std::string e_expr(expected_expr);
        std::string a_expr(actual_expr);

        failure_stream
            << setw(num_digits) << "i" << " "
            << setw(16) << (e_expr.size() <= 16 ? e_expr : "EXPECTED") << " "
            << setw(16) << (a_expr.size() <= 16 ? a_expr : "ACTUAL") << " "
            << setw(16) << "Difference" << "\n";

        // Loop through failed indices and print values
        for (Vec_Size_Type::const_iterator it = failed_indices.begin();
                it != end_indices;
                ++it)
        {
            if (std::fabs(expected[*it]) > 1.e-14)
                error = (actual[*it] - expected[*it]) / expected[*it];
            else
                error = actual[*it] - expected[*it];

            failure_stream
                << setw(num_digits) << *it << " "
                << setw(16) << expected[*it] << " "
                << setw(16) << actual[*it] << " "
                << setw(16) << error << "\n";
        }
        failure << failure_stream.str();

        return failure;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Vector equality comparison
 */
template<class T>
::testing::AssertionResult IsIterEq(
        const char* expected_expr,
        const char* actual_expr,
        const T* expected, const T* expected_end,
        const T* actual,   const T* actual_end)
{
    typedef std::size_t size_type;

    size_type expected_size = static_cast<size_type>(expected_end - expected);
    size_type actual_size   = static_cast<size_type>(actual_end - actual);

    // First, check that the sizes are equal
    if (expected_size != actual_size)
    {
        ::testing::AssertionResult failure = ::testing::AssertionFailure();

        failure << " Size of: " << actual_expr << "\n"
                << "  Actual: " << actual_size << "\n"
                << "Expected: " << expected_expr << ".size()\n"
                << "Which is: " << expected_size << "\n";
        return failure;
    }

    typedef std::vector<size_type> Vec_Size_Type;

    // Keep track of what elements failed
    Vec_Size_Type failed_indices;

    // Now loop through elements of the vectors and check for soft equivalence
    for (size_type i = 0; i < expected_size; ++i)
    {
        if (actual[i] != expected[i])
        {
            failed_indices.push_back(i);
        }
    }

    if (failed_indices.empty())
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        using std::setw;
        using std::setprecision;

        ::testing::AssertionResult failure = ::testing::AssertionFailure();

        failure << "Values in: " << actual_expr << "\n"
                << " Expected: " << expected_expr << "\n"
                << failed_indices.size() << " of "
                << expected_size << " elements differ:\n";

        // Only print the first 40 failures
        Vec_Size_Type::const_iterator end_indices = failed_indices.end();
        if (failed_indices.size() > 40)
        {
            failure << "(Truncating to first 40 failed values)\n";
            end_indices = failed_indices.begin() + 40;
        }

        // Calculate how many digits we need to space out
        unsigned int num_digits = calc_num_digits(failed_indices.back());

        // Construct our own stringstream because google test ignores setw
        std::ostringstream failure_stream;

        // Try to use user-given expressions for headers, but fall back if the
        // column length is exceeded
        std::string e_expr(expected_expr);
        std::string a_expr(actual_expr);

        failure_stream
            << setw(num_digits) << "i" << " "
            << setw(16) << (e_expr.size() <= 16 ? e_expr : "EXPECTED") << " "
            << setw(16) << (a_expr.size() <= 16 ? a_expr : "ACTUAL") << "\n";

        // Loop through failed indices and print values
        for (Vec_Size_Type::const_iterator it = failed_indices.begin();
                it != end_indices;
                ++it)
        {
            failure_stream
                << setw(num_digits) << *it << " "
                << setw(16) << expected[*it] << " "
                << setw(16) << actual[*it] << "\n";
        }
        failure << failure_stream.str();

        return failure;
    }
}

//---------------------------------------------------------------------------//
// >>> EXPLICIT INSTANTIATION
template ::testing::AssertionResult IsIterEq<int>(
        const char*, const char*,
        const int*, const int*,
        const int*, const int*);
template ::testing::AssertionResult IsIterEq<unsigned int>(
        const char*, const char*,
        const unsigned int*, const unsigned int*,
        const unsigned int*, const unsigned int*);
template ::testing::AssertionResult IsIterEq<float>(
        const char*, const char*,
        const float*, const float*,
        const float*, const float*);
template ::testing::AssertionResult IsIterEq<double>(
        const char*, const char*,
        const double*, const double*,
        const double*, const double*);

//---------------------------------------------------------------------------//
} // end namespace nemesis

//---------------------------------------------------------------------------//
//                 end of Gtest_Functions.cc
//---------------------------------------------------------------------------//
