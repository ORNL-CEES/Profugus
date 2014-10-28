//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/mc/Source_Transporter.hh
 * \author Thomas M. Evans
 * \date   Tue May 13 09:20:07 2014
 * \brief  Source_Transporter class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef core_mc_Source_Transporter_hh
#define core_mc_Source_Transporter_hh

#include <memory>

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "Source.hh"
#include "Domain_Transporter.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Source_Transporter
 * \brief Transport particles through geometry from source.
 *
 * This defines a Monte Carlo transporter that, given a source, transports
 * particles through phase-space and tallies results.  It serves as the
 * fundamental piece of either K-code or fixed-source MC solvers (a
 * fixed-source solver is generally a thin wrapper around this class).
 *
 * It solves the fixed source problem using a domain replication (DR) parallel
 * strategy.  In DR the entire mesh is replicated across all domains.
 */
/*!
 * \example mc/test/tstSource_Transporter.cc
 *
 * Test of Source_Transporter.
 */
//===========================================================================//

class Source_Transporter
{
  public:
    //@{
    //! Typedefs.
    typedef Domain_Transporter                   Transporter_t;
    typedef Source                               Source_t;
    typedef Transporter_t::Physics_t             Physics_t;
    typedef Transporter_t::Geometry_t            Geometry_t;
    typedef Transporter_t::SP_Physics            SP_Physics;
    typedef Transporter_t::SP_Geometry           SP_Geometry;
    typedef Transporter_t::SP_Particle           SP_Particle;
    typedef Transporter_t::SP_Variance_Reduction SP_Variance_Reduction;
    typedef Transporter_t::SP_Fission_Sites      SP_Fission_Sites;
    typedef Transporter_t::SP_Tallier            SP_Tallier;
    typedef std::shared_ptr<Source_t>            SP_Source;
    typedef Physics_t::RCP_Std_DB                RCP_Std_DB;
    typedef def::size_type                       size_type;
    //@}

  private:
    // >>> DATA

    // Problem geometry implementation.
    SP_Geometry d_geometry;

    // Problem physics implementation.
    SP_Physics d_physics;

    // Fixed source.
    SP_Source d_source;

    // Variance reduction.
    SP_Variance_Reduction d_var_reduction;

    // Tally controller.
    SP_Tallier d_tallier;

    // Domain transporter.
    Transporter_t d_transporter;

  public:
    // Constructor.
    Source_Transporter(RCP_Std_DB db, SP_Geometry geometry, SP_Physics physics);

    // Assign the source.
    void assign_source(SP_Source source);

    // Solve the fixed-source problem.
    void solve();

    // Set fission sampling.
    void sample_fission_sites(SP_Fission_Sites fis_sites, double keff);

    // Set the variance reduction.
    void set(SP_Variance_Reduction vr);

    // Set the tally controller
    void set(SP_Tallier tallier);

    // >>> ACCESSORS

    //! Get the tallies.
    SP_Tallier tallier() const { return d_tallier; }

    //! Get the transporter.
    const Transporter_t& transporter() const { return d_transporter; }

    //! Get the source.
    const Source_t& source() const { REQUIRE(d_source); return *d_source; }

  private:
    // >>> IMPLEMENTATION

    // Nodes and node id.
    int d_node, d_nodes;

    // Print out frequency for particle histories.
    double d_print_fraction;
    size_type d_print_count;
};

} // end namespace profugus

#endif // core_mc_Source_Transporter_hh

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.hh
//---------------------------------------------------------------------------//
