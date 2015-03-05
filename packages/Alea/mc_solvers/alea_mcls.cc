// Alea driver for MCLS solvers
#include "comm/global.hh"
#include "LinearSystem.hh"
#include "LinearSystemFactory.hh"
#include "AleaTypedefs.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include <MCLS_MCSASolverManager.hpp>
#include <MCLS_MultiSetLinearProblem.hpp>
#include <MCLS_TpetraAdapter.hpp>

using namespace alea;

int main( int argc, char *argv[] )
{
    // Initialize parallel communication.
    Teuchos::GlobalMPISession mpi_session( &argc, &argv );
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    Teuchos::FancyOStream out( Teuchos::rcpFromRef( std::cout ) );
    out.setOutputToRootOnly( 0 );
    out.setShowProcRank( true );

    // Timed execution block.
    {
	// Setup a total timer.
	Teuchos::RCP<Teuchos::Time> total_timer =
	    Teuchos::TimeMonitor::getNewCounter("Alea: Total");
	Teuchos::TimeMonitor total_monitor( *total_timer );

	// Read in command line options.
	std::string xml_input_filename;
	Teuchos::CommandLineProcessor clp(false);
	clp.setOption( "xml-in-file",
		       &xml_input_filename,
		       "The XML file to read into a parameter list" );
	clp.parse(argc,argv);

	// Build the parameter list from the xml input.
	Teuchos::RCP<Teuchos::ParameterList> plist =
	    Teuchos::rcp( new Teuchos::ParameterList() );
	Teuchos::updateParametersFromXmlFile(
	    xml_input_filename, Teuchos::inoutArg(*plist) );
	Teuchos::RCP<Teuchos::ParameterList> mcls_list = 
	    Teuchos::rcpFromRef( plist->sublist("MCLS",true) );

	// Build a communicator for the sets.
	int num_sets = mcls_list->get<int>("Number of Sets");
	int set_size = comm->getSize() / num_sets;
	int set_id = std::floor( Teuchos::as<double>(comm->getRank()) /
				 Teuchos::as<double>(set_size) );
	Teuchos::RCP<const Teuchos::Comm<int> > set_comm =
	    comm->split( set_id, comm->getRank() );

	// Get the raw set communicator.
	Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_set_comm = 
	    Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( set_comm );
	Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_set_comm = 
	    mpi_set_comm->getRawMpiComm();
	MPI_Comm raw_set_comm = (*opaque_set_comm)();
    
	// Set the default communicator with profugus as the set communicator.
	profugus::set_default( raw_set_comm );
	profugus::inherit( raw_set_comm );
    
	// Build the profugus linear system on each set.
	Teuchos::RCP<LinearSystem> system;
	{
	    Teuchos::RCP<Teuchos::Time> setup_timer =
		Teuchos::TimeMonitor::getNewCounter("SPN: Setup");
	    Teuchos::TimeMonitor setup_monitor( *setup_timer );
	    system =
		alea::LinearSystemFactory::buildLinearSystem(plist);
	}
	Teuchos::RCP<const MATRIX> pA = system->getMatrix();
	Teuchos::RCP<const MV> pb = system->getRhs();
	Teuchos::RCP<MV> px( new MV(pA->getDomainMap(),1) );

	// Extract the linear problem.
	Teuchos::RCP<const CRS_MATRIX> A =
	    Teuchos::rcp_dynamic_cast<const CRS_MATRIX>( pA );
	Teuchos::RCP<const VECTOR> b = pb->getVector( 0 );
	Teuchos::RCP<VECTOR> x = px->getVectorNonConst( 0 );
	Teuchos::RCP<MCLS::MultiSetLinearProblem<VECTOR,CRS_MATRIX> > problem =
	    Teuchos::rcp( 
		new MCLS::MultiSetLinearProblem<VECTOR,CRS_MATRIX>(
		    comm, num_sets, set_id, A, x, b) );

	// Build the MCLS solver.
	MCLS::MCSASolverManager<VECTOR,CRS_MATRIX> solver_manager( problem, mcls_list );

	// Solve the problem.
	solver_manager.solve();
    }
    
    // Output final timing.
    Teuchos::TableFormat& format = Teuchos::TimeMonitor::format();
    format.setPrecision(5);
    Teuchos::TimeMonitor::summarize();
    
    return 0;
}

