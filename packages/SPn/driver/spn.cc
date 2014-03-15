//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   driver/spn.cc
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
#include "comm/global.hh"
#include "utils/Definitions.hh"
#include "Manager.hh"

// Parallel specs.
int node  = 0;
int nodes = 0;

//---------------------------------------------------------------------------//
// Print instructions on how to run the neutronics executable

void print_usage()
{
    if (node == 0)
    {
        std::cout << "Usage: neutronics -i XMLFILE" << std::endl;
        std::cout << "Executes the neutronics executable using XMLFILE "
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

    // nodes
    node  = profugus::node();
    nodes = profugus::nodes();

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
        profugus::Manager manager;

        // setup the problem
        manager.setup(xml_file);
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

    profugus::finalize();
    return 0;
}

//---------------------------------------------------------------------------//
//                 end of spn.cc
//---------------------------------------------------------------------------//
