//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Arnoldi.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 14:41:30 2014
 * \brief  Arnoldi class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_Arnoldi_hh
#define solvers_Arnoldi_hh

#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "AnasaziTypes.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "EigenvalueSolver.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Arnoldi
 * \brief Compute dominant eigenvalue/eigenvector pair of operator.
 *
 * Given an Epetra operator, this class computes the eigenvalue of largest
 * magnitude and the corresponding eigenvector.  The Anasazi Block Krylov Schur
 * solver is currently used to solve the eigenproblem.
 *
 * \sa Arnoldi.cc for detailed descriptions.
 */
/*!
 * \example solvers/test/tstArnoldi.cc
 *
 * Test of Arnoldi.
 */
//===========================================================================//

class Arnoldi : public EigenvalueSolver<Epetra_MultiVector,Epetra_Operator>
{
  public:
    //@{
    //! Typedefs.
    typedef Epetra_MultiVector                            MV;
    typedef Epetra_Operator                               OP;
    typedef Anasazi::BasicEigenproblem<double,MV,OP>      Eigenproblem;
    typedef Anasazi::BlockKrylovSchurSolMgr<double,MV,OP> KrylovSchur;
    typedef Anasazi::Eigensolution<double,MV>             Eigensolution;

    typedef Teuchos::RCP<MV>           RCP_MV;
    typedef Teuchos::RCP<OP>           RCP_OP;
    typedef Teuchos::RCP<Eigenproblem> RCP_Eigenproblem;
    //@}

  private:

    // CMFD Operator
    RCP_OP d_A;

    // Teuchos Parameter List
    RCP_ParameterList d_pl;

  public:

    // Constructor
    Arnoldi( RCP_ParameterList db );

    // Set operator for Arnoldi
    void set_operator( RCP_OP op );

    // Set tolerance
    void set_tolerance(double tol)
    {
        Require(tol>0.0);
        b_tolerance = tol;
        d_pl->set("Convergence Tolerance",tol);
    }

    // Set iteration limit
    void set_max_iters(int iters)
    {
        Require(iters>0);
        b_max_iters = iters;
        d_pl->set("Maximum Restarts",iters);
    }

    // Solve eigenproblem
    void solve( double &eval, RCP_MV evec );
};

} // end namespace profugus

#endif // solvers_Arnoldi_hh

//---------------------------------------------------------------------------//
//                 end of Arnoldi.hh
//---------------------------------------------------------------------------//
