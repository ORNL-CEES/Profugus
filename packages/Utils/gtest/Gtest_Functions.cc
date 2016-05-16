//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/gtest/Gtest_Functions.cc
 * \author Seth R Johnson
 * \date   Tue Apr 02 11:41:41 2013
 * \brief  Gtest_Functions member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Gtest_Functions.hh"

#include <Utils/config.h>
#include <vector>
#include <iomanip>
#include <sstream>
#include <limits>
#include <cmath>
#include <cstddef>
#include <iterator>

#include "Utils/harness/Soft_Equivalence.hh"
#include "Utils/comm/Logger.hh"
#include "Utils/comm/global.hh"
#include "Utils/comm/Timer.hh"
#include "Utils/comm/Timing.hh"

#include "gtest.h"
#include "gtest-internals.hh"

namespace profugus
{
//---------------------------------------------------------------------------//
// INTERNAL-USE CLAASES
//---------------------------------------------------------------------------//
class ParallelHandler : public ::testing::EmptyTestEventListener
{
  public:
    // This will only be called after MPI_Init, so we can access these data
    ParallelHandler()
      : d_node(profugus::node())
      , d_num_nodes(profugus::nodes())
    {
        /* * */
    }

  public:
    virtual void OnTestProgramStart(const ::testing::UnitTest& /*unit_test*/)
    {
        // Write to cout
        profugus::logger().set("screen", &std::cout, profugus::DEBUG);

        profugus::log(profugus::DIAGNOSTIC)
            << "Testing on " << d_num_nodes << " processors";

#ifdef NEMESIS_PRINT_TEST_TIMING
        d_timer.start();
#endif
    }

#ifdef NEMESIS_PRINT_TEST_TIMING
    virtual void OnTestProgramEnd(const ::testing::UnitTest& /*unit_test*/)
    {
        // Process and output timing diagnostics
        profugus::global_barrier();
        d_timer.stop();
        double total_time = d_timer.TIMER_CLOCK();
        if (total_time > 0)
        {
            profugus::Timing_Diagnostics::report(std::cout, total_time);
        }

        // Output final timing
        if (d_node == 0)
        {
            std::cout << "Total execution time: "
                << std::scientific << std::setprecision(4) << total_time
                << " seconds.\n";
        }
    }
#endif

    virtual void OnTestStart(const ::testing::TestInfo& test_info)
    {
        // Barrier before starting the run
        profugus::global_barrier();
    }

    virtual void OnTestEnd(const ::testing::TestInfo& test_info)
    {
        // Barrier at the end of each test
        profugus::global_barrier();
    }

  private:
    int            d_node;
    int            d_num_nodes;
    profugus::Timer d_timer;
};

//---------------------------------------------------------------------------//
//! Class for printing failure messages on non-master nodes
class NonMasterResultPrinter : public ::testing::EmptyTestEventListener
{
    typedef ::testing::TestPartResult TestPartResult;

  public:
    // This will only be called after MPI_Init, so we can access these data
    NonMasterResultPrinter()
      : d_node(profugus::node())
    {
        /* * */
    }


    virtual void OnTestPartResult(const TestPartResult& result)
    {
        // If the test part succeeded, we don't need to do anything.
        if (result.type() == TestPartResult::kSuccess)
            return;

        ::testing::internal::ColoredPrintf(
                ::testing::internal::COLOR_RED, "[  FAILED  ] ");

        std::ostringstream os;

        if (result.file_name())
        {
            os << result.file_name() << ":";
        }
        if (result.line_number() >= 0)
        {
            os << result.line_number() << ":";
        }
        os << " Failure on node " << d_node << ":\n"
           << result.message();

        std::cout << os.str() << std::endl;
    }

  private:
    int d_node;
};


//---------------------------------------------------------------------------//
// NEMESIS FUNCTIONS
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

    ENSURE(num_digits > 0);
    return num_digits;
}

//-------------------------------------------------------------------------//
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
    profugus::initialize(argc, argv);
    const int node = profugus::node();

    // Initialize google test
    ::testing::InitGoogleTest(&argc, argv);

    // Gets hold of the event listener list.
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    if (node != 0)
    {
        // Don't print test completion messages (default pretty printer)
        delete listeners.Release(listeners.default_result_printer());
        // Instead, print just failure message (in case test fails on just one
        // node)
        listeners.Append(new NonMasterResultPrinter());
    }
    // Adds a listener to the end.  Google Test takes the ownership.
    listeners.Append(new ParallelHandler());

    // Run everything
    int failed = RUN_ALL_TESTS();

    // Accumulate the result so that all processors will have the same result
    profugus::global_sum(failed);

    // If no tests were run, there's a problem.
    if (testing::UnitTest::GetInstance()->test_to_run_count() == 0)
    {
        if (node == 0)
        {
            ::testing::internal::ColoredPrintf(
                    ::testing::internal::COLOR_RED,
                    "[  FAILED  ] ");
            std::cout << "No tests are written/enabled!" << std::endl;
        }

        failed = 1;
    }

    // Finalize MPI
    profugus::global_barrier();
    profugus::finalize();

    // Print final results
    if (node == 0)
    {
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
} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Gtest_Functions.cc
//---------------------------------------------------------------------------//
