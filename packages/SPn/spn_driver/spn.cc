//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn_driver/spn.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 12 22:24:55 2014
 * \brief  SPn Mini-App executable.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <string>
#include <iostream>
#include <algorithm>

#include "harness/DBC.hh"
#include "comm/Timer.hh"
#include "comm/global.hh"
#include "comm/Timing_Diagnostics.hh"
#include "comm/P_Stream.hh"
#include "utils/Definitions.hh"
#include "Manager.hh"
#include "SPn/config.h"

#include <Teuchos_TimeMonitor.hpp>

// Parallel specs.
int node  = 0;
int nodes = 0;

//---------------------------------------------------------------------------//
// Print instructions on how to run the neutronics executable

void print_usage()
{
    if (node == 0)
    {
        std::cout << "Usage: xspn -i XMLFILE" << std::endl;
        std::cout << "Executes the xspn executable using XMLFILE "
                  << "as the input file." << std::endl;
        exit(1);
    }
}

//---------------------------------------------------------------------------//
// Parse the input arguments

std::string parse_input_arguments(const def::Vec_String &arguments)
{
    // First, search for "-h" or "--help"
    if (std::find(arguments.begin(), arguments.end(), "-h") != arguments.end()
        || std::find(arguments.begin(), arguments.end(), "--help") !=
       arguments.end())
    {
        print_usage();
    }

    // Search for "-i"
    auto iter = std::find(arguments.begin(), arguments.end(), "-i");
    if (iter == arguments.end() || iter == arguments.end()-1)
    {
        if (node == 0)
        {
            std::cout << std::endl << "ERROR: Missing xml input filename."
                      << std::endl << std::endl;
        }
        print_usage();
    }

    // Get the xml filename
    const std::string &xml_filename = *(iter+1);
    if (xml_filename.empty())
    {
        if (node == 0)
        {
            std::cout << std::endl << "ERROR: Missing xml input filename."
                      << std::endl << std::endl;
        }
        print_usage();
    }

    return xml_filename;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    profugus::initialize(argc, argv);

    // start timing
    profugus::global_barrier();
    profugus::Timer timer;
    timer.start();

    // nodes
    node  = profugus::node();
    nodes = profugus::nodes();

    profugus::pcout << "=======================================\n"
                    << "    Profugus SPN Mini-APP              \n"
                    << "    (C) ORNL, Battelle, 2014           \n"
                    << "=======================================\n"
                    << profugus::endl;

    // process input arguments
    def::Vec_String arguments(argc - 1);
    std::string     xml_file;
    for (int c = 1; c < argc; c++)
    {
        arguments[c - 1] = argv[c];
    }
    xml_file = parse_input_arguments(arguments);

    try
    {
        // make the manager
        spn::Manager manager;

        // setup the problem
        manager.setup(xml_file);

        // solve the problem
        manager.solve();

        // output
        manager.output();
    }
    catch (const profugus::assertion &a)
    {
        std::cout << "Caught profugus assertion " << a.what() << std::endl;
        exit(1);
    }
    catch (const std::exception &a)
    {
        std::cout << "Caught standard assertion " << a.what() << std::endl;
        exit(1);
    }
    catch (...)
    {
        std::cout << "Caught assertion of unknown origin." << std::endl;
        exit(1);
    }

    // process and output timing diagnostics
    profugus::global_barrier();
    timer.stop();
    double total = timer.TIMER_CLOCK();
    profugus::Timing_Diagnostics::report(std::cout, total);

    // output final timing
    profugus::pcout << "\n" << "Total execution time : "
                    << profugus::scientific
                    << profugus::setprecision(4)
                    << total << " seconds." << profugus::endl << profugus::endl;

#ifdef USE_TRILINOS_TIMING
    // output final timing from trilinos components
    Teuchos::TableFormat &format = Teuchos::TimeMonitor::format();
    format.setPrecision(5);
    Teuchos::TimeMonitor::summarize();
#endif

    profugus::finalize();
    return 0;
}

//---------------------------------------------------------------------------//
//                 end of spn.cc
//---------------------------------------------------------------------------//
