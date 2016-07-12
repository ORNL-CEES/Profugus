//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/KCode_Solver.cuh
 * \author Steven Hamilton
 * \date   Mon May 19 10:30:32 2014
 * \brief  KCode_Solver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_KCode_Solver_cuh
#define cuda_mc_KCode_Solver_cuh

#include <thrust/device_vector.h>

#include "cuda_utils/CudaDBC.hh"
#include "Definitions.hh"
#include "Solver.hh"
#include "Fission_Source.cuh"
#include "Tallier.cuh"
#include "Keff_Tally.cuh"
#include "Source_Transporter.hh"

namespace cuda_mc 
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

template <class Geometry>
class KCode_Solver : public Solver<Geometry>
{
    typedef Solver<Geometry> Base;

  public:
    //@{
    //! Typedefs.
    typedef Geometry                                Geometry_t;
    typedef Source_Transporter<Geometry_t>          Source_Transporter_t;
    typedef Teuchos::RCP<Teuchos::ParameterList>    RCP_Std_DB;
    typedef std::shared_ptr<Source_Transporter_t>   SP_Source_Transporter;
    typedef Fission_Source<Geometry_t>              FS_t;
    typedef Tallier<Geometry_t>                     Tallier_t;
    typedef Keff_Tally<Geometry_t>                  Keff_Tally_t;
    typedef std::vector<Fission_Site>               Host_Fission_Sites;
    typedef thrust::device_vector<Fission_Site>     Dev_Fission_Sites;
    typedef std::shared_ptr<FS_t>                   SP_Fission_Source;
    typedef std::shared_ptr<Host_Fission_Sites>     SP_Host_Fission_Sites;
    typedef std::shared_ptr<Dev_Fission_Sites>      SP_Dev_Fission_Sites;
    typedef std::shared_ptr<Tallier_t>              SP_Tallier;
    typedef cuda::Shared_Device_Ptr<Tallier_t>      SDP_Tallier;
    typedef std::shared_ptr<Keff_Tally_t>           SP_Keff_Tally;
    typedef cuda::Shared_Device_Ptr<Keff_Tally_t>   SDP_Keff_Tally;
    typedef def::size_type                          size_type;

  private:

    // >>> DATA

    SDP_Keff_Tally d_keff_tally;
    SP_Keff_Tally d_keff_tally_host;

    // Problem database.
    RCP_Std_DB d_db;

    // Source transporter.
    SP_Source_Transporter d_transporter;

    // Fission source.
    SP_Fission_Source       d_source;
    SP_Host_Fission_Sites   d_host_sites;
    SP_Dev_Fission_Sites    d_dev_sites;

    // Inactive tallier
    SDP_Tallier d_inactive_tallier;

  public:

    // Constructor.
    KCode_Solver(RCP_Std_DB db);

    // Set the underlying fixed-source transporter and fission source.
    void set(SP_Source_Transporter transporter,
             SP_Fission_Source     source,
             SP_Tallier            tallier);

    // >>> ACCESSORS

    //! Number of (inactive or active) cycles run so far
    unsigned int num_cycles() const
    {
        return d_keff_tally_host->cycle_count();
    }

    // Get keff value
    double keff() const { return d_keff_tally_host->mean(); }

    /*
     * \brief Access inactive tallier
     *
     * This is provided (to be used before solve()) so that extra tallies
     * beside k-effective can be added to the inactive cycles.
     */
    SDP_Tallier inactive_tallier()
    {
        REQUIRE(d_build_phase >= ASSIGNED);
        REQUIRE(d_inactive_tallier.get_host_ptr());
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

    using Base::b_tallier;

    // Build phase for the solver
    Build_Phase d_build_phase;

    // Whether to suppress printing (true for non-master node)
    bool d_quiet;

    // Number of particles per cycle (constant weight).
    double d_Np;

    // Number of particles in batch
    size_type d_batch_size;
};

} // end namespace cuda_mc

#endif // cuda_mc_KCode_Solver_cuh

//---------------------------------------------------------------------------//
//                 end of KCode_Solver.cuh
//---------------------------------------------------------------------------//

