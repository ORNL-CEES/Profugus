//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Eigenvalue_Solver.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Mon Mar 10 14:20:33 2014
 * \brief  Eigenvalue_Solver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Eigenvalue_Solver_hh
#define spn_Eigenvalue_Solver_hh

#include <SPn/config.h>

#include "Ifpack_Preconditioner.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"

// ML has to be optional for Windows compatibility
#ifdef USE_ML
#include "ml_MultiLevelPreconditioner.h"
#endif

#include "comm/Timer.hh"
#include "solvers/EigenvalueSolver.hh"
#include "Solver_Base.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Eigenvalue_Solver
 * \brief Solve the SPN equations for an eigenvalue problem.
 *
 * This class solves the SPN equations for the following eigenvalue problem:
 * \f[
   \mathbf{A}\mathbf{u} = \frac{1}{k}\mathbf{B}\mathbf{u}\:,
 * \f]
 * where \b A is the SPN operator, \f$\mathbf{u}=\{u_1,u_2,u_3,u_4\}\f$ are
 * the transformed moments of the flux, \b B is the fission matrix, and \e k
 * is the (dominant) eigenvalue.  In order to write the scalar flux (0\e th
 * SPN moment) into the state use write_state().
 *
 * \sa spn::Linear_System
 */
/*!
 * \example spn/test/tstEigenvalue_Solver.cc
 *
 * Test of Eigenvalue_Solver.
 */
//===========================================================================//

class Eigenvalue_Solver : public Solver_Base
{
  public:
    //@{
    //! Typedefs.
    typedef Solver_Base                        Base;
    typedef Epetra_MultiVector                 MV;
    typedef Epetra_Operator                    OP;
    typedef Teuchos::RCP<OP>                   RCP_Epetra_Op;
    typedef profugus::EigenvalueSolver<EPETRA> Eigensolver;
    typedef Teuchos::RCP<Eigensolver>          RCP_Eigensolver;
    //@}

  private:
    // >>> DATA

    // Eigenvector.
    RCP_Vector d_u;

    // Eigenvalue.
    double d_keff;

    // Eigensolver.
    RCP_Eigensolver d_eigensolver;

    // Preconditioners.
    // We have to hold on to these because they are given to
    // other classes as raw pointers and we need to make sure
    // they don't get destoyed (until they should)
    Teuchos::RCP<Ifpack_Preconditioner> d_ifpack_prec;
#ifdef USE_ML
    Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> d_ml_prec;
#endif

  public:
    // Constructor.
    explicit Eigenvalue_Solver(RCP_ParameterList db);

    // Set up the solver.
    void setup(RCP_Dimensions dim, RCP_Mat_DB mat, RCP_Mesh mesh,
               RCP_Indexer indexer, RCP_Global_Data data);

    // Solve the SPN eigenvalue equations.
    void solve();

    // Write the scalar-flux (eigenvector) into the state.
    void write_state(State_t &state);

    // >>> ACCESSORS

    //! Get eigenvalue (\e keff).
    double get_eigenvalue() const { return d_keff; }

    //! Get eigen-vector (in transformed \e u space).
    const Vector_t& get_eigenvector() const { return *d_u; }

  private:
    // >>> IMPLEMENTATION

    // Set db defaults
    void set_default_parameters();

    // Build the preconditioner.
    RCP_Epetra_Op build_preconditioner(RCP_Dimensions dim, RCP_Mat_DB mat,
                                       RCP_Mesh mesh, RCP_Indexer indexer,
                                       RCP_Global_Data data);

    // Timer.
    profugus::Timer d_timer;
};

} // end namespace profugus

#endif // spn_Eigenvalue_Solver_hh

//---------------------------------------------------------------------------//
//                 end of Eigenvalue_Solver.hh
//---------------------------------------------------------------------------//
