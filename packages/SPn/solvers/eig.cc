//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/eig.cc
 * \author Steveb Hamilton
 * \date   Wed Mar 12 22:24:55 2014
 * \brief  Executable to solve eigenvalue problem
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <string>
#include <iostream>
#include <algorithm>

#include "harness/DBC.hh"
#include "comm/Timer.hh"
#include "comm/global.hh"
#include "comm/P_Stream.hh"
#include "utils/Definitions.hh"
#include "SPn/config.h"

#include <Teuchos_TimeMonitor.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "EpetraExt_CrsMatrixIn.h"

#include "LinAlgTypedefs.hh"
#include "PreconditionerBuilder.hh"
#include "EigenvalueSolverBuilder.hh"


// Parallel specs.
int node  = 0;
int nodes = 0;

//---------------------------------------------------------------------------//
// Print instructions on how to run the neutronics executable

void print_usage()
{
    if (node == 0)
    {
        std::cout << "Usage: xeig -i XMLFILE" << std::endl;
        std::cout << "Executes the xeig executable using XMLFILE "
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

void run_problem(const std::string &xml_file)
{
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // make the master parameterlist
    auto pl = Teuchos::rcp(new Teuchos::ParameterList(""));

    // read the data on every domain
    Teuchos::updateParametersFromXmlFileAndBroadcast(
        xml_file.c_str(), pl.ptr(), *comm);

    auto node = Teuchos::rcp( new profugus::TpetraTypes::NODE() );

    typedef profugus::EpetraTypes   MatType;
    typedef MatType::MATRIX         MATRIX;

#ifdef COMM_MPI
    Epetra_MpiComm ecomm(profugus::communicator);
#else
    Epetra_SerialComm ecomm;
#endif

    // Get name of A matrix and read
    std::cout << "Loading A from file." << std::endl;
    auto a_file = pl->get<std::string>("A_filename");
    MATRIX *A_ptr;
    int err;
    err = EpetraExt::MatrixMarketFileToCrsMatrix(a_file.c_str(),ecomm,A_ptr);
    INSIST(err == 0, "Error Reading A");
    Teuchos::RCP<MATRIX> A(A_ptr);

    // Get name of B matrix and read
    std::cout << "Loading B from file." << std::endl;
    auto b_file = pl->get<std::string>("B_filename");
    MATRIX *B_ptr;
    err = EpetraExt::MatrixMarketFileToCrsMatrix(b_file.c_str(),ecomm,B_ptr);
    INSIST(err == 0, "Error Reading B");
    Teuchos::RCP<MATRIX> B(B_ptr);

    // Build preconditioner
    std::cout << "Building preconditioner." << std::endl;
    profugus::Timer timer;
    timer.start();
    auto P = profugus::PreconditionerBuilder<MatType>::build_preconditioner(
        A, pl);

    timer.stop();
    std::cout << "Preconditioner construction took " << timer.TIMER_CLOCK()
        << " seconds." << std::endl;

    // Build solver
    auto solver = profugus::EigenvalueSolverBuilder<MatType>::build_solver(
        pl, A, B, P);


    // Build Eigenvector
    auto x = Teuchos::rcp( new MatType::MV( A->RowMap(), 1 ) );

    timer.reset();
    timer.start();

    // Solve
    double lambda = 1000.0;
    std::cout << "Solving." << std::endl;
    solver->solve( lambda, x );

    timer.stop();
    double total = timer.TIMER_CLOCK();

    std::cout << "Computed eigenvalue of " << lambda <<
        " in " << total << " seconds" << std::endl;
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

    run_problem(xml_file);

    // process and output timing diagnostics
    profugus::global_barrier();
    timer.stop();
    double total = timer.TIMER_CLOCK();

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
