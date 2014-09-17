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

template <class T>
class Eigenvalue_Solver : public Solver_Base_Tmpl<T>
{
  public:
    //@{
    //! Typedefs.
    typedef Solver_Base_Tmpl<T>                         Base;
    typedef typename T::MV                              MV;
    typedef typename T::OP                              OP;
    typedef Teuchos::RCP<OP>                            RCP_OP;
    typedef profugus::EigenvalueSolver<T>               Eigensolver;
    typedef Teuchos::RCP<Eigensolver>                   RCP_Eigensolver;
    typedef Linear_System<T>                            Linear_System_t;
    typedef typename Linear_System_t::External_Source   External_Source;
    typedef typename Linear_System_t::RCP_Timestep      RCP_Timestep;
    typedef typename Linear_System_t::Vector_t          Vector_t;
    typedef typename Linear_System_t::RCP_ParameterList RCP_ParameterList;
    typedef typename Linear_System_t::RCP_Vector        RCP_Vector;
    typedef typename Linear_System_t::RCP_Dimensions    RCP_Dimensions;
    typedef typename Linear_System_t::RCP_Mat_DB        RCP_Mat_DB;
    typedef typename Linear_System_t::RCP_Mesh          RCP_Mesh;
    typedef typename Linear_System_t::RCP_Indexer       RCP_Indexer;
    typedef typename Linear_System_t::RCP_Global_Data   RCP_Global_Data;
    //@}

    using Base::b_db;
    using Base::b_system;

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
    void write_state(State &state);

    // >>> ACCESSORS

    //! Get eigenvalue (\e keff).
    double get_eigenvalue() const { return d_keff; }

    //! Get eigen-vector (in transformed \e u space).
    Teuchos::RCP<const Vector_t> get_eigenvector() const { return d_u; }

  private:
    // >>> IMPLEMENTATION

    // Set db defaults
    void set_default_parameters();

    // Build the preconditioner.
    RCP_OP build_preconditioner(RCP_Dimensions dim, RCP_Mat_DB mat,
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
