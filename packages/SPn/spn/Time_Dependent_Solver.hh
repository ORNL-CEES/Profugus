//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Time_Dependent_Solver.hh
 * \author Thomas M. Evans
 * \date   Fri Apr 04 00:09:50 2014
 * \brief  Time_Dependent_Solver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Time_Dependent_Solver_hh
#define spn_Time_Dependent_Solver_hh

#include "comm/Timer.hh"
#include "solvers/StratimikosSolver.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "Solver_Base.hh"
#include "VectorTraits.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Time_Dependent_Solver
 * \brief Solve a time-dependent SPN problem.
 *
 * This is just a stub for future time-dependent development.
 */
/*!
 * \example spn/test/tstTime_Dependent_Solver.cc
 *
 * Test of Time_Dependent_Solver.
 */
//===========================================================================//

template <class T>
class Time_Dependent_Solver : public Solver_Base_Tmpl<T>
{
  public:
    //@{
    //! Typedefs.
    typedef Solver_Base_Tmpl<T>                         Base;
    typedef typename T::MV                              MV;
    typedef typename T::OP                              OP;
    typedef StratimikosSolver<T>                        Linear_Solver_t;
    typedef Linear_System<T>                            Linear_System_t;
    typedef typename Linear_System_t::External_Source   External_Source;
    typedef typename Linear_System_t::RCP_Timestep      RCP_Timestep;
    typedef typename Linear_System_t::RCP_ParameterList RCP_ParameterList;
    typedef typename Linear_System_t::RCP_MV            RCP_MV;
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

    // Solver.
    Linear_Solver_t d_solver;

    // Solution vector.
    RCP_MV d_lhs;

    // Timestep controller.
    RCP_Timestep d_dt;

  public:
    // Constructor.
    explicit Time_Dependent_Solver(RCP_ParameterList db);

    // Set up the solver.
    void setup(RCP_Dimensions dim, RCP_Mat_DB mat, RCP_Mesh mesh,
               RCP_Indexer indexer, RCP_Global_Data data, bool adjoint = false);

    // Solve the SPN equations.
    void solve(Teuchos::RCP<const External_Source> q);

    // Write the scalar-flux into the state.
    void write_state(State &state);

    // Write problem to file
    void write_problem_to_file() const;

    // >>> ACCESSORS

    //! Get LHS solution vector (in transformed \e u space).
    Teuchos::RCP<const MV> get_LHS() const { return d_lhs; }

  private:
    // >>> IMPLEMENTATION

    // Timer.
    profugus::Timer d_timer;
};

} // end namespace profugus

#endif // spn_Time_Dependent_Solver_hh

//---------------------------------------------------------------------------//
//                 end of Time_Dependent_Solver.hh
//---------------------------------------------------------------------------//
