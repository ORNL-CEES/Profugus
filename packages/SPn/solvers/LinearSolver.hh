//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/LinearSolver.hh
 * \author Steven Hamilton
 * \date   Mon Jul 01 11:48:09 2013
 * \brief  LinearSolver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_LinearSolver_hh
#define solvers_LinearSolver_hh

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "harness/DBC.hh"
#include "harness/Warnings.hh"
#include "comm/P_Stream.hh"
#include "utils/String_Functions.hh"

#include "LinAlgTypedefs.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class LinearSolver
 * \brief Base class for solving linear system with templated
 *  vector/operator interface.
 */
//===========================================================================//

template <LinAlgType T>
class LinearSolver
{
  public:
    //@{
    //! Typedefs.
    typedef typename LinAlgTypedefs<T>::MV  MV;
    typedef typename LinAlgTypedefs<T>::OP  OP;
    typedef Teuchos::ParameterList          ParameterList;
    typedef Teuchos::RCP<ParameterList>     RCP_ParameterList;
    //@}

  protected:

    enum Verbosity { NONE=0, LOW=1, MEDIUM=2, HIGH=3, DEBUG=4 };

  public:

    // Constructor
    LinearSolver(RCP_ParameterList db)
        : b_db(db)
    {
        // Get stopping tolerance off of DB or set a default.
        b_tolerance = db->get<double>("tolerance", 1.0e-6);

        // Get max iterations off of db or set a default.
        b_max_iters = db->get<int>("max_itr", 100);

        // Get verbosity or set default
        b_verbosity = LOW;
        if( db->isParameter("verbosity") )
        {
            std::string verb =
                profugus::to_lower(db->get<std::string>("verbosity"));
            if( verb=="none" )
            {
                b_verbosity = NONE;
            }
            else if( verb=="low" )
            {
                b_verbosity = LOW;
            }
            else if( verb=="medium" )
            {
                b_verbosity = MEDIUM;
            }
            else if( verb=="high" )
            {
                b_verbosity = HIGH;
            }
            else if( verb=="debug" )
            {
                b_verbosity = DEBUG;
            }
        }
    }

    // Virtual destructor
    virtual ~LinearSolver(){};

    // Set operator
    virtual void set_operator( Teuchos::RCP<OP> A )
    {
        REQUIRE( !A.is_null() );
        b_A = A;
    }

    // Set preconditioner
    virtual void set_preconditioner( Teuchos::RCP<OP> P )
    {
        ADD_WARNING("Preconditioning not supported by " << b_label);
    }

    // Solve
    virtual void solve( Teuchos::RCP<MV>       x,
                        Teuchos::RCP<const MV> b ) = 0;

    // Set tolerance
    virtual void set_tolerance(double tol)
    {
        REQUIRE(tol>0.0);
        b_tolerance = tol;
    }

    // Set max iterations
    virtual void set_max_iters(int iters)
    {
        REQUIRE(iters>0);
        b_max_iters=iters;
    }

    // Return iteration count from latest solve
    virtual int num_iters() const { return b_num_iters; }

    // Did last solver meet convergence tolerance?
    virtual bool converged() const { return b_converged; }

    // Return solver label
    virtual const std::string & solver_label() const { return b_label; }

    // Print status to screen
    virtual void print_status( Verbosity level, int iter, double res )
    {
        if( b_verbosity >= level )
        {
            profugus::pout << b_label << " at iteration " << iter
                << " has residual norm " << res << profugus::endl;
        }
    }

  protected:

    RCP_ParameterList b_db;
    Teuchos::RCP<OP>  b_A;
    double            b_tolerance;
    int               b_max_iters;
    int               b_num_iters;
    bool              b_converged;
    std::string       b_label;
    Verbosity         b_verbosity;

};

} // end namespace profugus

#endif // solvers_LinearSolver_hh

//---------------------------------------------------------------------------//
//              end of solvers/LinearSolver.hh
//---------------------------------------------------------------------------//
