//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Source_Transporter.hh
 * \author Stuart Slattery
 * \date   Tue May 13 09:20:07 2014
 * \brief  Source_Transporter class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Source_Transporter_hh
#define cuda_mc_Source_Transporter_hh

#include <memory>

#include "cuda_utils/Shared_Device_Ptr.hh"
#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "Source.hh"
#include "Domain_Transporter.hh"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Source_Transporter
 * \brief Transport particles through geometry from source.
 *
 * This defines a Monte Carlo transporter that, given a source, transports
 * particles through phase-space and tallies results.  It serves as the
 * fundamental piece of either K-code or fixed-source CUDA_MC solvers (a
 * fixed-source solver is generally a thin wrapper around this class).
 *
 * It solves the fixed source problem using a domain replication (DR) parallel
 * strategy.  In DR the entire mesh is replicated across all domains.
 */
/*!
 * \example cuda_mc/test/tstSource_Transporter.cc
 *
 * Test of Source_Transporter.
 */
//===========================================================================//

template <class Geometry>
class Source_Transporter
{
  public:
    //@{
    //! Typedefs.
    typedef Domain_Transporter<Geometry>                  Transporter_t;
    typedef Source<Geometry>                              Source_t;
    typedef typename Transporter_t::Physics_t             Physics_t;
    typedef typename Transporter_t::Geometry_t            Geometry_t;
    typedef typename Transporter_t::SDP_Physics           SDP_Physics;
    typedef typename Transporter_t::SDP_Geometry          SDP_Geometry;
    typedef typename Transporter_t::SDP_Particle_Vector   SDP_Particle_Vector;
    typedef typename Transporter_t::SP_Fission_Sites      SP_Fission_Sites;
    typedef typename Transporter_t::SP_Tallier            SP_Tallier;
    typedef typename Transporter_t::SP_Variance_Reduction SP_Variance_Reduction;
    typedef std::shared_ptr<Source_t>                     SP_Source;
    typedef Teuchos::RCP<Teuchos::ParameterList>          RCP_Std_DB;
    typedef def::size_type                                size_type;
    //@}

  private:
    // >>> DATA

    // Problem geometry implementation.
    SDP_Geometry d_geometry;

    // Problem physics implementation.
    SDP_Physics d_physics;

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
    Source_Transporter(const RCP_Std_DB& db, 
                       const SDP_Geometry& geometry, 
                       const SDP_Physics& physics);

    // Assign the source.
    void assign_source(const SP_Source& source);

    // Solve the fixed-source problem.
    void solve();

    // Set fission sampling.
    void sample_fission_sites(const SP_Fission_Sites& fis_sites, double keff);

    // Set the variance reduction.
    void set(SP_Variance_Reduction vr);

    // Set the tally controller
    void set(const SP_Tallier& tallier);

    // >>> ACCESSORS

    //! Get the tallies.
    SP_Tallier tallier() const { return d_tallier; }

    //! Get the transporter.
    const Transporter_t& transporter() const { return d_transporter; }

    //! Get the source.
    const Source_t& source() const { REQUIRE(d_source); return *d_source; }

    //! Get the geometry.
    SDP_Geometry geometry() const { REQUIRE(d_geometry); return d_geometry; }

    //! Get the physics.
    SDP_Physics physics() const { REQUIRE(d_physics); return d_physics; }

  private:
    // >>> IMPLEMENTATION

    // Nodes and node id.
    int d_node, d_nodes;

    // Print out frequency for particle histories.
    double d_print_fraction;
    size_type d_print_count;

    // particle vector size.
    size_type d_vector_size;

    // number of batches.
    int d_num_batch;
};

} // end namespace cuda_profugus

#endif // cuda_mc_Source_Transporter_hh

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.hh
//---------------------------------------------------------------------------//
