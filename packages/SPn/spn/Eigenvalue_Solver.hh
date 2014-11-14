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

#include "harness/DBC.hh"
#include "comm/Timer.hh"
#include "solvers/EigenvalueSolver.hh"
#include "solvers/LinAlgTypedefs.hh"
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
    typedef Teuchos::RCP<MV>                            RCP_MV;
    typedef profugus::EigenvalueSolver<T>               Eigensolver;
    typedef Teuchos::RCP<Eigensolver>                   RCP_Eigensolver;
    typedef Linear_System<T>                            Linear_System_t;
    typedef Teuchos::RCP<Linear_System_t>               RCP_Linear_System;
    typedef typename Linear_System_t::External_Source   External_Source;
    typedef typename Linear_System_t::RCP_Timestep      RCP_Timestep;
    typedef typename Linear_System_t::RCP_ParameterList RCP_ParameterList;
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
    RCP_MV d_u;

    // Eigenvalue.
    double d_keff;

    // Eigensolver.
    RCP_Eigensolver d_eigensolver;

    // Material database
    RCP_Mat_DB d_mat;

  public:
    // Constructor.
    explicit Eigenvalue_Solver(RCP_ParameterList db);

    // Set up the solver.
    void setup(RCP_Dimensions dim, RCP_Mat_DB mat, RCP_Mesh mesh,
               RCP_Indexer indexer, RCP_Global_Data data, bool adjoint = false);

    // Set up the solver from a pre-built linear system.
    void setup(RCP_Mat_DB mat, RCP_Mesh mesh, RCP_Indexer indexer,
               RCP_Global_Data data, RCP_Linear_System system,
               bool adjoint = false);

    // Solve the SPN eigenvalue equations.
    void solve(Teuchos::RCP<const External_Source> q);

    // Write the scalar-flux (eigenvector) into the state.
    void write_state(State &state);

    // >>> ACCESSORS

    //! Get eigenvalue (\e keff).
    double get_eigenvalue() const { return d_keff; }

    //! Get eigen-vector (in transformed \e u space).
    Teuchos::RCP<const MV> get_eigenvector() const { return d_u; }

    //! Write problem matrices to file
    void write_problem_to_file() const;

  private:
    // >>> IMPLEMENTATION

    // Set db defaults
    void set_default_parameters();

    // Build the preconditioner.
    RCP_OP build_preconditioner(RCP_Dimensions dim, RCP_Mat_DB mat,
                                RCP_Mesh mesh, RCP_Indexer indexer,
                                RCP_Global_Data data);

    // Apply adjoint (transpose) to the operators.
    void apply_transpose(bool adjoint)
    {
        NOT_IMPLEMENTED(
            "Failed to apply transpose on non-Epetra implementations.");
    }

    // Timer.
    profugus::Timer d_timer;
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Apply transpose on Epetra operators.
 */
template<>
inline void Eigenvalue_Solver<EpetraTypes>::apply_transpose(bool adjoint)
{
    // set transpose for adjoint calculations
    b_system->get_Operator()->SetUseTranspose(adjoint);
    b_system->get_fission_matrix()->SetUseTranspose(adjoint);
}

} // end namespace profugus

#endif // spn_Eigenvalue_Solver_hh

//---------------------------------------------------------------------------//
//                 end of Eigenvalue_Solver.hh
//---------------------------------------------------------------------------//
