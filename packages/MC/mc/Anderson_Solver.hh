//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Anderson_Solver.hh
 * \author Thomas M. Evans
 * \date   Tue Apr 07 20:32:28 2015
 * \brief  Anderson_Solver class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Anderson_Solver_hh
#define MC_mc_Anderson_Solver_hh

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "comm/global.hh"
#include "solvers/AndersonSolver.hh"
#include "Solver.hh"
#include "Anderson_Operator.hh"
#include "Keff_Tally.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Anderson_Solver
 * \brief Solve the eigenvalue problem using Anderson.
 */
/*!
 * \example mc/test/tstAnderson_Solver.cc
 *
 * Test of Anderson_Solver.
 */
//===========================================================================//

template<class T>
class Anderson_Solver : public Solver
{
  public:
    typedef Anderson_Operator<T>                     Operator;
    typedef Teuchos::RCP<Teuchos::ParameterList>     RCP_Std_DB;
    typedef typename Operator::RCP_MAP               RCP_MAP;
    typedef Teuchos::RCP<Operator>                   RCP_Operator;
    typedef AndersonSolver<T>                        Anderson_t;
    typedef std::shared_ptr<Anderson_t>              SP_Anderson;
    typedef typename Operator::SP_Source_Transporter SP_Source_Transporter;
    typedef typename Operator::SP_Fission_Source     SP_Fission_Source;
    typedef std::shared_ptr<Keff_Tally>              SP_Keff_Tally;

  private:
    // >>> DATA

    // Problem database.
    RCP_Std_DB d_db;

    // Anderson operator.
    RCP_Operator d_operator;

    // Keff tally.
    SP_Keff_Tally d_keff_tally;

    // Anderson non-linear solver.
    SP_Anderson d_anderson;

  public:
    // Constructor.
    Anderson_Solver(RCP_Std_DB db);

    // Set the underlying fixed-source transporter and fission source.
    void set(SP_Source_Transporter transporter, SP_Fission_Source source);

    // >>> INHERITED INTERFACE

    // Solve the problem.
    void solve() override;

    // Call to reset the solver and tallies for another calculation.
    void reset() override;

    // >>> ACCESSORS

    //! Keff tally for tracking k-effective.
    SP_Keff_Tally keff_tally() const { return d_keff_tally; }

  private:
    // >>> IMPLEMENTATION

    // Initialize the problem.
    void initialize();

    // Run inactive cycle(s).
    void run_inactive();

    // Do Anderson solve.
    void anderson_solve();

    // Run active cycles.
    void run_active();

    // Number of particles per cycle (constant weight).
    double d_Np;

    // Node and nodes
    int d_nodes, d_node;

    // Anderson k-eff value.
    double d_k_anderson;

    // Set-constant communicator.
    profugus::Communicator_t d_set_comm;
};

} // end namespace profugus

#endif // MC_mc_Anderson_Solver_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Anderson_Solver.hh
//---------------------------------------------------------------------------//
