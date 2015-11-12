//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MultiSplitting.cc
 * \author Massimiliano Lupo Pasini
 * \brief  Perform Multi Splitting with either deterministic Richardson 
 * iterations or MC updates.
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "Teuchos_ArrayRCP.hpp"

#include "MultiSplitting.hh"
#include "LinearSolverFactory.hh"
#include "harness/DBC.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param A  Problem matrix
 * \param pl ParameterList
 *
 * Behavior is controlled by following PL entries on the nested
 * "MultiSplitting"
 * sublist:
 *  - block_size(int)             : >0 (1)
 *  - overlapping(SCALAR)         : >0.0 (0.0)
 */
//---------------------------------------------------------------------------//

MultiSplitting::MultiSplitting( Teuchos::RCP<Teuchos::ParameterList> &pl )
{
    d_pl = pl;        
             
    Teuchos::RCP<NODE> node = KokkosClassic::Details::getNode<NODE>();
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();            

    Teuchos::RCP<Teuchos::ParameterList> prob_pl =
        Teuchos::sublist(pl,"Problem");

    std::string filename = prob_pl->get<std::string>("matrix_filename");
    d_A = Tpetra::MatrixMarket::Reader<CRS_MATRIX>::readSparseFile(
        filename,comm,node);
	
    Teuchos::RCP<const MAP> map = d_A->getDomainMap();
        
    std::string rhs_file = prob_pl->get<std::string>("rhs_filename");
    d_b = Tpetra::MatrixMarket::Reader<CRS_MATRIX>::readVectorFile(
        rhs_file,comm,node,map);
    
    std::cout<<"Is upper triangular: "<<d_A->isUpperTriangular()<<std::endl;
    std::cout<<"Is lower triangular: "<<d_A->isLowerTriangular()<<std::endl;

    d_multisplitting = Teuchos::RCP<LinearSystem_MultiSplitting>( new LinearSystem_MultiSplitting(d_pl, d_A, d_b));
      
    // Get MultiSplitting sublist
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(d_pl,"MultiSplitting");             
             
    d_inner_solver = d_multisplitting->getInnerSolverType();

    d_num_blocks = d_multisplitting->getNumBlocks(); 

    VALIDATE( d_inner_solver == "richardson" || d_inner_solver == "monte_carlo",
         "the only iterative solvers admitted are Richardson and MonteCarlo" );                  
    d_divergence_tol = mat_pl->get("divergence_tolerance",1.0e4);
    
    b_max_iterations = mat_pl->get<int>("max_iterations");

    b_tolerance = mat_pl->get<MAGNITUDE>("tolerance");
    
    if( mat_pl->isType<std::string>("verbosity") )
    {
        std::string verbosity = mat_pl->get<std::string>("verbosity", "low");
        VALIDATE(verbosity=="none"   || verbosity=="low"  ||
                 verbosity=="medium" || verbosity=="high" ||
                 verbosity=="debug",
                 "Invalid verbosity specified.");
        if( verbosity == "none")
        {
            b_verbosity = NONE;
        }
        else if( verbosity == "low" )
        {
            b_verbosity = LOW;
        }
        else if( verbosity == "medium" )
        {
            b_verbosity = MEDIUM;
        }
        else if( verbosity == "high" )
        {
            b_verbosity = HIGH;
        }
        else if( verbosity == "debug" )
        {
            b_verbosity = DEBUG;
        }
    }
     
}


//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

Teuchos::RCP<MV> MultiSplitting::computeIteration()
{
    Teuchos::RCP<const MAP> map = d_A->getDomainMap();
    Teuchos::RCP<MV> x_partial = Teuchos::RCP<MV>( new MV(map, d_num_blocks) );
    Teuchos::RCP<MV> x_average;

    for(unsigned int p=0; p!=10; ++p)
    {
        std::cout<<"p= "<<p<<std::endl;
        splitting split= d_multisplitting->buildSplitting(p);
        Teuchos::RCP<CRS_MATRIX> A = split.A;
        Teuchos::RCP<MV> b = split.b;
        Teuchos::RCP<MV> E = split.E;
        Teuchos::RCP<MV> x_p;
        
//    	Teuchos::RCP<alea::AleaSolver> solver =
//       		alea::LinearSolverFactory::buildSolver(d_inner_solver,A,d_pl);
//    	solver->apply(*b,*x_p);
    }
    std::cout<<"Esco del ciclo for"<<std::endl;
    return x_average;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Perform MultiSplitting.
 */
//---------------------------------------------------------------------------//
void MultiSplitting::solve(Teuchos::RCP<MV> &x) 
{

    // For now we only support operating on a single vector
    REQUIRE( x->getNumVectors() == 1 );
    // Compute initial residual
    MV r( d_b->getMap(),1 );
    r.update(1.0,*(d_b),0.0);

    Teuchos::ArrayRCP<MAGNITUDE> r_norm(1);
    r.norm2(r_norm());
    if( r_norm[0] == 0.0 )
        return;
        
    MAGNITUDE r0_norm = r_norm[0];

    b_num_iters = 0;
    while( true )
    {
        b_num_iters++;

        // Compute residual r = b - A*x
        d_A->apply(*x,r);
        r.update(1.0,*d_b,-1.0);
        std::cout<<"residual norm: "<<r_norm[0]<<std::endl;
        Teuchos::RCP<const MV> r_pointer = Teuchos::RCP<const MV>(new MV(r));

        // update the multisplitting system with the residual as RHS
        d_multisplitting.reset( new LinearSystem_MultiSplitting(d_pl, d_A, r_pointer) );

	//std::cout<<d_inner_solver<<std::endl;
        Teuchos::RCP<MV> aux = computeIteration();
        // Check convergence on true (rather than preconditioned) residual
        r.norm2(r_norm());
 
	std::cout<<"verbosity: "<<b_verbosity<<std::endl;
        if( b_verbosity >= HIGH )
        {
            std::cout << "Relative residual norm at iteration " << b_num_iters
                << " is " << r_norm[0]/r0_norm << std::endl;
        }

        // Check for convergence
        if( r_norm[0]/r0_norm < b_tolerance )
        {
            if( b_verbosity >= LOW )
            {
                std::cout << "Richardson Iteration converged after "
                    << b_num_iters << " iterations." << std::endl;
            }
            break;
        }

        // Check for max iterations
        if( b_num_iters >= b_max_iterations )
        {
            std::cout << "Richardson Iteration reached maximum iteration "
                 << "count with relative residual norm of "
                 << r_norm[0]/r0_norm << std::endl;
            b_num_iters = -1;
            break;
        }

        // Check for divergence
        if( r_norm[0]/r0_norm > d_divergence_tol)
        {
            if( b_verbosity >= LOW )
            {
                std::cout << "Richardson Iteration diverged after "
                    << b_num_iters << " iterations." << std::endl;
            }
            b_num_iters = -b_num_iters;
            break;
        }

    }
         
}

} // namespace alea

