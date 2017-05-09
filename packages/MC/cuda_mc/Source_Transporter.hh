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
#include "Particle_Vector.cuh"

namespace cuda_mc
{

// Forward declarations to avoid including Cuda headers
template <class Geom> class Physics;
template <class Geom> class Domain_Transporter_DMM;
template <class Geom> class Tallier_DMM;
template <class Geom> class VR_Roulette;
template <class Geom> class Source;

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
    typedef Physics<Geometry>                             Physics_t;
    typedef Domain_Transporter_DMM<Geometry>              Transporter_DMM_t;
    typedef Tallier_DMM<Geometry>                         Tallier_DMM_t;
    typedef VR_Roulette<Geometry>                         VR_Roulette_t;
    typedef cuda::Shared_Device_Ptr<Geometry>             SDP_Geometry;
    typedef cuda::Shared_Device_Ptr<Physics_t>            SDP_Physics;
    typedef cuda::Shared_Device_Ptr<VR_Roulette_t>        SDP_VR;
    typedef std::shared_ptr<Transporter_DMM_t>            SP_Transporter_DMM;
    typedef std::shared_ptr<Tallier_DMM_t>                SP_Tallier_DMM;
    typedef Teuchos::RCP<Teuchos::ParameterList>          RCP_Std_DB;
    typedef def::size_type                                size_type;
    typedef thrust::device_vector<Fission_Site>           Fission_Site_Vector;
    typedef std::shared_ptr<Fission_Site_Vector>          SP_Fission_Site_Vec;
    typedef Source<Geometry>                              Source_t;
    typedef std::shared_ptr<Source_t>                     SP_Source;
    typedef Particle_Vector_DMM<Geometry>                 Particle_Vector_DMM_t;
    typedef std::shared_ptr<Particle_Vector_DMM_t>        SP_Particle_Vec_DMM;
    //@}

  private:
    // >>> DATA

    // Problem geometry implementation.
    SDP_Geometry d_geometry;

    // Problem physics implementation.
    SDP_Physics d_physics;

    // Tally controller.
    SP_Tallier_DMM d_tallier;

    // Variance reduction.
    SDP_VR d_vr;

    // Domain transporter.
    SP_Transporter_DMM d_transporter;

    // Particle_Vector
    SP_Particle_Vec_DMM d_particle_vec;

  public:
    // Constructor.
    Source_Transporter(RCP_Std_DB   db,
                       SDP_Geometry geometry,
                       SDP_Physics  physics);

    ~Source_Transporter()
    {
        if (d_verbosity >= LOW)
        {
            std::cout << "ST spent " << d_source_time
                << " generating source particles, " << d_transport_time
                << " transporting particles and " << d_sort_time
                << " sorting particles" << std::endl;
        }
    }

    // Set tallier
    void set(SP_Tallier_DMM tallier);

    // Solve the fixed-source problem.
    void solve(SP_Source source) const;

    // Set fission sampling.
    void sample_fission_sites(SP_Fission_Site_Vec fis_sites, double keff);

    //! Get number of fission sites created during last transport.
    int num_sampled_fission_sites();

    // >>> ACCESSORS

    //! Get the tallies.
    SP_Tallier_DMM tallier() const { return d_tallier; }

    //! Get the geometry.
    SDP_Geometry geometry() const { return d_geometry; }

    //! Get the physics.
    SDP_Physics physics() const { return d_physics; }

    //! Get the transporter.
    SP_Transporter_DMM transporter() const { return d_transporter; }

  private:
    // >>> IMPLEMENTATION

    enum Sort_Type {ALIVE, MATID, GROUP, CELL};
    enum Verbosity {NONE, LOW, MEDIUM, HIGH};

    // Nodes and node id.
    int d_node, d_nodes;
    int d_block_size;

    Sort_Type d_sort_type;
    Verbosity d_verbosity;

    mutable double d_source_time, d_sort_time, d_transport_time;
};

} // end namespace cuda_mc

#endif // cuda_mc_Source_Transporter_hh

//---------------------------------------------------------------------------//
//                 end of Source_Transporter.hh
//---------------------------------------------------------------------------//
