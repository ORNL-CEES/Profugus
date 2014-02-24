//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Davidson_Eigensolver.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 14:41:41 2014
 * \brief  Davidson_Eigensolver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_Davidson_Eigensolver_hh
#define solvers_Davidson_Eigensolver_hh

#include "AnasaziGeneralizedDavidsonSolMgr.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "AnasaziEpetraAdapter.hpp"

#include "EigenvalueSolver.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Davidson_Eigensolver
 * \brief Solve k-eigenvalue problem using Generalized Davidson solver
 */
/*!
 * \example solvers/test/tstDavidson_Eigensolver.cc
 *
 * Test of Davidson_Eigensolver.
 */
//===========================================================================//

class Davidson_Eigensolver :
        public EigenvalueSolver<Epetra_MultiVector,Epetra_Operator>
{
    typedef EigenvalueSolver<Epetra_MultiVector,Epetra_Operator> Base;

  public:
    //@{
    //! Typedefs.
    typedef Epetra_MultiVector                 MV;
    typedef Epetra_Operator                    OP;
    typedef Teuchos::RCP<OP>                   RCP_OP;
    typedef Anasazi::MultiVecTraits<double,MV> MultiVecTraits;
    //@}

  public:
    //! Constructor
    Davidson_Eigensolver( RCP_ParameterList db,
                          RCP_OP            LHS,
                          RCP_OP            RHS );

    //! Register preconditioner with solver
    void set_preconditioner( RCP_OP prec )
    {
        Require( prec != Teuchos::null );
        d_prec = prec;
    }

    //! Perform setup operations
    void setup();

    //! Solve eigenproblem
    void solve( double           &lambda,
                Teuchos::RCP<MV>  x );

  private:

    // Solver database
    RCP_ParameterList d_db;

    // Operators
    RCP_OP d_LHS;
    RCP_OP d_RHS;
    RCP_OP d_prec;
};

} // end namespace profugus

#endif // solvers_Davidson_Eigensolver_hh

//---------------------------------------------------------------------------//
//                 end of Davidson_Eigensolver.hh
//---------------------------------------------------------------------------//
