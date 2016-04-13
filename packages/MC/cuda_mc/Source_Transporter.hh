//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Source_Transporter.hh
 * \author Thomas M. Evans
 * \date   Tue May 13 09:20:07 2014
 * \brief  Source_Transporter class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Source_Transporter_hh
#define cuda_mc_Source_Transporter_hh

#include <memory>
#include <thrust/device_vector.h>

#include "utils/Definitions.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Shared_Device_Ptr can be included in host-only code
// No "raw" cuda headers should appear here
#include "cuda_utils/Shared_Device_Ptr.hh"
#include "cuda_utils/CudaDBC.hh"

#include "Definitions.hh"

namespace cuda_mc
{

// Forward declarations to avoid including Cuda headers
template <class Geom> class Particle;
template <class Geom> class Physics;
template <class Geom> class Domain_Transporter;
template <class Geom> class Tallier;
template <class Geom> class VR_Roulette;
class RNG_Control;

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

template <class Geometry>
class Source_Transporter
{
  public:
    //@{
    //! Typedefs.
    typedef Particle<Geometry>                            Particle_t;
    typedef Physics<Geometry>                             Physics_t;
    typedef Domain_Transporter<Geometry>                  Transporter_t;
    typedef Tallier<Geometry>                             Tallier_t;
    typedef VR_Roulette<Geometry>                         VR_Roulette_t;
    typedef cuda::Shared_Device_Ptr<Geometry>             SDP_Geometry;
    typedef cuda::Shared_Device_Ptr<Particle_t>           SDP_Particle;
    typedef cuda::Shared_Device_Ptr<Physics_t>            SDP_Physics;
    typedef cuda::Shared_Device_Ptr<Transporter_t>        SDP_Transporter;
    typedef cuda::Shared_Device_Ptr<Tallier_t>            SDP_Tallier;
    typedef cuda::Shared_Device_Ptr<VR_Roulette_t>        SDP_VR;
    typedef std::shared_ptr<RNG_Control>                  SP_RNG_Control;
    typedef Teuchos::RCP<Teuchos::ParameterList>          RCP_Std_DB;
    typedef def::size_type                                size_type;
    typedef thrust::device_vector<Fission_Site>           Fission_Site_Vector;
    typedef std::shared_ptr<Fission_Site_Vector>          SP_Fission_Site_Vec;
    //@}

  private:
    // >>> DATA

    // Problem geometry implementation.
    SDP_Geometry d_geometry;

    // Problem physics implementation.
    SDP_Physics d_physics;

    // Tally controller.
    SDP_Tallier d_tallier;

    // Variance reduction.
    SDP_VR d_vr;

    // Domain transporter.
    SDP_Transporter d_transporter;

    // RNG Control
    SP_RNG_Control d_rng_control;

  public:
    // Constructor.
    Source_Transporter(RCP_Std_DB   db,
                       SDP_Geometry geometry,
                       SDP_Physics  physics,
                       SDP_Tallier  tallier = SDP_Tallier());

    // Solve the fixed-source problem.
    template <class Src_Type>
    void solve(std::shared_ptr<Src_Type> source) const;

    // Set fission sampling.
    void sample_fission_sites(SP_Fission_Site_Vec fis_sites, double keff);

    // >>> ACCESSORS

    //! Get the tallies.
    SDP_Tallier tallier() const { return d_tallier; }

    //! Get the geometry.
    SDP_Geometry geometry() const { return d_geometry; }

    //! Get the physics.
    SDP_Physics physics() const { return d_physics; }

    //! Get the transporter.
    SDP_Transporter transporter() const { return d_transporter; }

  private:
    // >>> IMPLEMENTATION

    // Nodes and node id.
    int d_node, d_nodes;
};

} // end namespace cuda_mc

#endif // cuda_mc_Source_Transporter_hh

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.hh
//---------------------------------------------------------------------------//
