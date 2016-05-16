//----------------------------------*-C++-*----------------------------------// /*!
/*!
 * \file   SPn/spn/Fixed_Source_Solver.hh
 * \author Thomas M. Evans
 * \date   Mon Feb 17 21:00:33 2014
 * \brief  Fixed_Source_Solver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_spn_Fixed_Source_Solver_hh
#define SPn_spn_Fixed_Source_Solver_hh

#include "comm/Timer.hh"
#include "solvers/StratimikosSolver.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "Solver_Base.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Fixed_Source_Solver
 * \brief Solve the SPN equations for a fixed source problem.
 *
 * This class solves the SPN equations for the following fixed-source problem:
 * \f[
   \mathbf{A}\mathbf{u} = \mathbf{Q}\:,
 * \f]
 * where \b A is the SPN operator, \f$\mathbf{u}=\{u_1,u_2,u_3,u_4\}\f$ are
 * the transformed moments of the flux, and \b Q is the source.  The system is
 * solved for an input external source.  In order to write the scalar flux
 * (0\e th SPN moment) into the state use write_state().
 *
 * \sa spn::Linear_System,  profugus::Isotropic_Source
 */
/*!
 * \example spn/test/tstFixed_Source_Solver.cc
 *
 * Test of Fixed_Source_Solver.
 */
//===========================================================================//

template <class T>
class Fixed_Source_Solver : public Solver_Base_Tmpl<T>
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

  public:
    // Constructor.
    explicit Fixed_Source_Solver(RCP_ParameterList db);

    // Set up the solver.
    void setup(RCP_Dimensions dim, RCP_Mat_DB mat, RCP_Mesh mesh,
               RCP_Indexer indexer, RCP_Global_Data data, bool adjoint = false);

    // Solve the SPN equations.
    void solve(Teuchos::RCP<const External_Source> q);

    // Write the scalar-flux into the state.
    void write_state(State &state);

    // Write problem matrix to file
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

#endif // SPn_spn_Fixed_Source_Solver_hh

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.hh
//---------------------------------------------------------------------------//
