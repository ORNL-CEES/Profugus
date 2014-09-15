//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/Fixed_Source_Solver.hh
 * \author Thomas M. Evans
 * \date   Mon Feb 17 21:00:33 2014
 * \brief  Fixed_Source_Solver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_tpetra_Fixed_Source_Solver_hh
#define spn_tpetra_Fixed_Source_Solver_hh

#include "comm/Timer.hh"
#include "solvers/StratimikosSolver.hh"
#include "Solver_Base.hh"

#include "solvers/LinAlgTypedefs.hh"

namespace profugus
{
namespace tpetra
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
 * \example spn_tpetra/test/tstFixed_Source_Solver.cc
 *
 * Test of Fixed_Source_Solver.
 */
//===========================================================================//

class Fixed_Source_Solver : public Solver_Base
{
  public:
    //@{
    //! Typedefs.
    typedef Solver_Base                              Base;
    typedef Linear_System_t::External_Source         External_Source;
    typedef Linear_System_t::RCP_Timestep            RCP_Timestep;
    typedef LinAlgTypedefs<TPETRA>            LinAlgImpl;
    typedef typename LinAlgImpl::OP           OP;
    typedef typename LinAlgImpl::MV           MV;
    typedef StratimikosSolver<TPETRA>                Linear_Solver_t;
    //@}

  private:
    // >>> DATA

    // Solver.
    Linear_Solver_t d_solver;

    // Solution vector.
    RCP_Vector d_lhs;

  public:
    // Constructor.
    explicit Fixed_Source_Solver(RCP_ParameterList db);

    // Set up the solver.
    void setup(RCP_Dimensions dim, RCP_Mat_DB mat, RCP_Mesh mesh,
               RCP_Indexer indexer, RCP_Global_Data data);

    // Solve the SPN equations.
    void solve(const External_Source &q);

    // Write the scalar-flux into the state.
    void write_state(State_t &state);

    // >>> ACCESSORS

    //! Get LHS solution vector (in transformed \e u space).
    Teuchos::RCP<const Vector_t> get_LHS() const { return d_lhs; }

  private:
    // >>> IMPLEMENTATION

    // Timer.
    profugus::Timer d_timer;
};

} // end namespace profugus
} // end namespace tpetra

#endif // spn_tpetra_Fixed_Source_Solver_hh

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.hh
//---------------------------------------------------------------------------//
