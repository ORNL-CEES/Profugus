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

#include "Keff_Solver.hh"

#include "comm/global.hh"
#include "solvers/AndersonSolver.hh"
#include "Anderson_Operator.hh"

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
class Anderson_Solver : public Keff_Solver
{
    typedef Keff_Solver Base;

  public:
    typedef Anderson_Operator<T>                     Operator;
    typedef Teuchos::RCP<Teuchos::ParameterList>     RCP_Std_DB;
    typedef typename Operator::RCP_MAP               RCP_MAP;
    typedef Teuchos::RCP<Operator>                   RCP_Operator;
    typedef AndersonSolver<T>                        Anderson_t;
    typedef std::shared_ptr<Anderson_t>              SP_Anderson;
    typedef typename Operator::SP_Source_Transporter SP_Source_Transporter;
    typedef typename Base::SP_Fission_Source         SP_Fission_Source;

  private:
    // >>> DATA

    // Problem database.
    RCP_Std_DB d_db;

    // Anderson operator.
    RCP_Operator d_operator;

    // Anderson non-linear solver.
    SP_Anderson d_anderson;

    // Bring base class types in.
    using Base::b_keff_tally;

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

    //! Get acceleration.
    SP_FM_Acceleration acceleration() const
    {
        return typename Base::SP_FM_Acceleration();
    }

    // >>> ACCESSORS

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
