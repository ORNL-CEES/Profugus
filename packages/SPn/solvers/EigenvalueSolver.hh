//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/EigenvalueSolver.hh
 * \author Thomas M. Evans, Steven P. Hamilton
 * \date   Fri Feb 21 14:33:00 2014
 * \brief  EigenvalueSolver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_EigenvalueSolver_hh
#define solvers_EigenvalueSolver_hh

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "harness/DBC.hh"
#include "utils/String_Functions.hh"

#include "LinAlgTypedefs.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class EigenvalueSolver
 * \brief Base class for solving linear system with templated
 *  vector/operator interface.
 */
//===========================================================================//

template <class T>
class EigenvalueSolver
{
  public:
    //@{
    //! Typedefs.
    typedef typename T::MV                  MV;
    typedef typename T::OP                  OP;
    typedef Teuchos::ParameterList          ParameterList;
    typedef Teuchos::RCP<ParameterList>     RCP_ParameterList;
    //@}

  protected:

    enum Verbosity { NONE=0, LOW=1, MEDIUM=2, HIGH=3, DEBUG=4 };

  public:

    // Constructor
    EigenvalueSolver(RCP_ParameterList db)
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
    virtual ~EigenvalueSolver(){};

    // Set operator
    virtual void set_operator( Teuchos::RCP<OP> A )
    {
        REQUIRE( !A.is_null() );
        b_A = A;
    }

    // Solve
    virtual void solve( double           &lambda,
                        Teuchos::RCP<MV>  x ) = 0;

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

#endif // solvers_EigenvalueSolver_hh

//---------------------------------------------------------------------------//
//                 end of EigenvalueSolver.hh
//---------------------------------------------------------------------------//
