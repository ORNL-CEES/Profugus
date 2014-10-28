//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/mc/KCode_Solver.hh
 * \author Thomas M. Evans
 * \date   Mon May 19 10:30:32 2014
 * \brief  KCode_Solver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef core_mc_KCode_Solver_hh
#define core_mc_KCode_Solver_hh

#include "Solver.hh"

#include "harness/DBC.hh"
#include "Source_Transporter.hh"
#include "Keff_Tally.hh"
#include "Fission_Source.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class KCode_Solver
 * \brief Solve an eigenvalue problem using the MC K-Code method.
 *
 * The total number of cycles in the k-code calculation is broken down as
 * follows:
 * \f[
   \mbox{num\_cycles} = \mbox{num\_active\_cycles} +
 * \mbox{num\_inactive\_cycles}\:.
 * \f]
 *
 * We obtain our tally controller from the source transporter, but we
 * temporarily disable all those tallies during the inactive cycles.
 *
 * \section kcode_solver_db Standard DB Entries for KCode_Solver
 *
 * The database entries that control the profugus::KCode_Solver are:
 *
 * \arg \c keff_init (double) initial keff value to use for fission source
 * initialization (default: 1.0)
 *
 * \arg \c num_cycles (int) total number of k-code cycles (default: 50)
 *
 * \arg \c num_inactive_cycles (int) number of inactive cycles to run before
 * accumulating statistics (default: 10)
 */
/*!
 * \example mc/test/tstKCode_Solver.cc
 *
 * Test of KCode_Solver.
 */
//===========================================================================//

class KCode_Solver : public Solver
{
  public:
    //@{
    //! Typedefs.
    typedef Source_Transporter                    Source_Transporter_t;
    typedef Source_Transporter_t::RCP_Std_DB      RCP_Std_DB;
    typedef std::shared_ptr<Source_Transporter_t> SP_Source_Transporter;
    typedef std::shared_ptr<Keff_Tally>           SP_Keff_Tally;
    typedef std::shared_ptr<Fission_Source>       SP_Fission_Source;
    typedef Fission_Source::SP_Fission_Sites      SP_Fission_Sites;

  private:
    // >>> DATA

    // Problem database.
    RCP_Std_DB d_db;

    // Source transporter.
    SP_Source_Transporter d_transporter;

    // Fission source.
    SP_Fission_Source d_source;
    SP_Fission_Sites  d_fission_sites;

    // Eigenvalue tally
    SP_Keff_Tally d_keff_tally;

    // Inactive tallier
    SP_Tallier d_inactive_tallier;

  public:
    // Constructor.
    KCode_Solver(RCP_Std_DB db);

    // Set the underlying fixed-source transporter and fission source.
    void set(SP_Source_Transporter transporter, SP_Fission_Source source);

    // >>> ACCESSORS

    //! Keff tally for tracking k-effective (read-only please!)
    SP_Keff_Tally keff_tally() const { return d_keff_tally; }

    //! Number of (inactive or active) cycles run so far
    auto num_cycles() const -> decltype(d_keff_tally->cycle_count())
    {
        return d_keff_tally->cycle_count();
    }

    /*
     * \brief Access inactive tallier
     *
     * This is provided (to be used before solve()) so that extra tallies
     * beside k-effective can be added to the inactive cycles.
     */
    SP_Tallier inactive_tallier()
    {
        REQUIRE(d_build_phase >= ASSIGNED);
        ENSURE(d_inactive_tallier);
        return d_inactive_tallier;
    }

    // >>> INHERITED INTERFACE

    // Solve the kcode problem.
    void solve();

    // Call to reset the solver and tallies for another kcode run
    void reset();

    // >>> PUBLIC INTERFACE

    // Call at the beginning of a solve
    void initialize();

    // Call for each inactive/active cycle
    void iterate();

    // Call to transition from inactive to active cycles
    void begin_active_cycles();

    // Call to finalize tallies
    void finalize();

  private:
    // >>> IMPLEMENTATION

    // Phases of construction, for error checking
    enum Build_Phase {
        CONSTRUCTED = 0,//!< after construction is complete
        ASSIGNED,       //!< after assigning transporter, source
        INACTIVE_SOLVE, //!< after the call to initialize()
        ACTIVE_SOLVE,   //!< after the call to begin_active_cycles()
        FINALIZED       //!< after the call to finalize()
    };

    // Build phase for the solver
    Build_Phase d_build_phase;

    // Whether to suppress printing (true for non-master node)
    bool d_quiet;

    // Number of particles per cycle (constant weight).
    double d_Np;
};

} // end namespace profugus

#endif // core_mc_KCode_Solver_hh

//---------------------------------------------------------------------------//
//                 end of KCode_Solver.hh
//---------------------------------------------------------------------------//
