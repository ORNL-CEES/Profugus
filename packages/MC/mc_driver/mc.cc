//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_driver/mc.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 12 22:24:55 2014
 * \brief  SPn Mini-App executable.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <hpx/hpx_init.hpp>
#include <boost/program_options.hpp>

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

//---------------------------------------------------------------------------//
int hpx_main( boost::program_options::variables_map& vm )
{
    profugus::pcout << "=======================================\n"
                    << "    Profugus MC Mini-APP               \n"
                    << "    (C) ORNL, Battelle, 2014           \n"
                    << "=======================================\n"
                    << profugus::endl;

    // start timing
    profugus::global_barrier();
    profugus::Timer timer;
    timer.start();

    // Get the XML input file.
    std::string xml_file = vm["input"].as<std::string>();

    // Run the problem.
    try
    {
        // make the manager
        mc::Manager manager;

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
                    << total << " seconds." << profugus::endl;

    // Finalize hpx.
    return hpx::finalize();
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[])
{
    // initialize mpi
    profugus::initialize(argc, argv);

    // process input arguments
    boost::program_options::options_description
	desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");
    desc_commandline.add_options()
	( "input",
	  boost::program_options::value<std::string>()->default_value("mc.xml"),
	  "Profugus XML input file");

    // initalize hpx. will call hpx_main and run the problem.
    int hpx_result = hpx::init( desc_commandline, argc, argv );

    // finalize mpi
    profugus::finalize();

    // Return the result from hpx.
    return hpx_result;
}

//---------------------------------------------------------------------------//
//                 end of mc.cc
//---------------------------------------------------------------------------//
