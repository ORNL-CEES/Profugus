//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fission_Source.cuh
 * \author Steven Hamilton
 * \date   Mon May 05 14:22:46 2014
 * \brief  Fission_Source class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Fission_Source_cuh
#define cuda_mc_Fission_Source_cuh

#include <thrust/device_vector.h>

#include "cuda_geometry/Cartesian_Mesh.hh"
#include "Fission_Rebalance.hh"
#include "Physics.cuh"
#include "Particle.cuh"
#include "RNG_Control.cuh"
#include "Source.cuh"
#include "Definitions.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Fission_Source
 * \brief Defines both an initial and within-cycle fission source for k-code
 * calculations.
 *
 * Through the database passed into the fission source the client specifies
 * the number of fission particles to run each cycle.  On cycles other than
 * the first, the number of fission source particles is determined by the
 * number of sites sampled during the previous generation.  To keep the
 * effective number of particles constant between cycles, fission particles
 * are born with the following weight,
 * \f[
   w = N_p / M\:,
 * \f]
 * where \f$N_p\f$ is the number of requested particles per cycle and \f$M\f$
 * is the total number of particles sampled from the fission source.  In the
 * first cycle \f$M = N_p\f$; however, in various parallel decompositions the
 * number of particles run in the first cycle may be slightly less than the
 * number of particles requested (but never more).  Thus, even in the first
 * cycle the weight of particles may be slightly greater than one because the
 * weight per particle is normalized to the \b requested number of particles
 * per cycle.  Using this weight preserves the constant total weight per
 * cycle.  The requested number of particles can be queried through the Np()
 * member function.
 *
 * \section fission_source_db Standard DB Entries for Fission_Source
 *
 * The database entries that control the profugus::Fission_Source are:
 *
 * \arg \c init_fission_src (vector<double>) box defined as (lox, hix, loy,
 * hiy, loz, hiz) in which the initial fission source is sampled (default: the
 * entire geometry)
 *
 * \arg \c Np (int) number of particles to use in each cycle (default:
 * 1000)
 */
/*!
 * \example mc/test/tstFission_Source.cc
 *
 * Test of Fission_Source.
 */
//===========================================================================//

template <class Geometry>
class Fission_Source : public Source<Geometry>
{
  public:
    //@{
    //! Typedefs.
    typedef Geometry                                    Geometry_t;
    typedef Physics<Geometry_t>                         Physics_t;
    typedef thrust::device_vector<Fission_Site>         Fission_Site_Vector;
    typedef Particle<Geometry_t>                        Particle_t;
    typedef cuda_utils::Space_Vector                    Space_Vector;
    typedef RNG_Control::RNG_State_t                    RNG_State_t;
    typedef cuda::Shared_Device_Ptr<Geometry>           SDP_Geometry;
    typedef cuda::Shared_Device_Ptr<Physics_t>          SDP_Physics;
    typedef std::shared_ptr<Fission_Rebalance>          SP_Fission_Rebalance;
    typedef Teuchos::RCP<Teuchos::ParameterList>        RCP_Std_DB;
    typedef def::size_type                              size_type;
    typedef std::shared_ptr<Fission_Site_Vector>        SP_Fission_Site_Vec;
    //@}

  private:
    // >>> DATA

    // Fission site container.
    SP_Fission_Site_Vec d_fission_site_vec;
    Fission_Site       *d_fission_sites;

    // Fission rebalance (across sets).
    SP_Fission_Rebalance d_fission_rebalance;

    // Physics
    SDP_Physics d_physics_host;
    Physics_t  *d_physics;

  public:
    // Constructor.
    Fission_Source(RCP_Std_DB db, SDP_Geometry geometry, SDP_Physics physics);

    // Build the initial fission source.
    void build_initial_source();

    // Build a source from a fission site container.
    void build_source(SP_Fission_Site_Vec &fission_sites);

    // Is this the initial source
    bool is_initial_source() const { return !d_have_sites; }

    // >>> DERIVED PUBLIC INTERFACE

    // Get a particle from the source.
    __device__ inline Particle_t get_particle(
        std::size_t idx, RNG_State_t *rng) const;

    //! Number of particles to transport in the source on the current domain.
    __host__ __device__
    size_type num_to_transport() const { return d_np_domain; }

    //! Total number of particles to transport in the entire problem/cycle.
    __host__ __device__ 
    size_type total_num_to_transport() const { return d_np_total; }

    // >>> CLASS ACCESSORS

    //! Get the current fission site container.
    SP_Fission_Site_Vec fission_sites() const
    {
        return d_fission_site_vec;
    }

    //! Total number of requested particles per cycle.
    __host__ __device__ size_type Np() const { return d_np_requested; }

    //! Set a new number per cycle.
    void update_Np(size_type np) { d_np_requested = np; }

  private:
    // >>> IMPLEMENTATION

    typedef Source<Geometry> Base;
    using Base::b_geometry;

    // Build the domain replicated fission source.
    void build_DR();

    // Sample the geometry.
    __device__ inline
    int sample_geometry(Space_Vector &r, const Space_Vector &omega,
                        Particle_t &p, RNG_State_t *rng) const;

    int d_nodes;

    // Have we populated fission_sites yet?
    bool d_have_sites;

    // Initial fission source lower coords and width.
    Space_Vector d_lower;
    Space_Vector d_width;

    // Requested particles per cycle.
    size_type d_np_requested;

    // Number of particles: total, domain
    size_type d_np_total, d_np_domain;

    // Particle weight.
    double d_wt;
};

} // end namespace cuda_mc

#include "Fission_Source.i.cuh"

#endif // cuda_mc_Fission_Source_cuh

//---------------------------------------------------------------------------//
//                 end of Fission_Source.cukh
//---------------------------------------------------------------------------//
