//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Domain_Transporter.hh
 * \author Stuart Slattery
 * \date   Mon May 12 12:02:13 2014
 * \brief  Domain_Transporter class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Domain_Transporter_hh
#define cuda_mc_Domain_Transporter_hh

#include <memory>

#include "Physics.hh"
#include "Variance_Reduction.hh"
#include "Tallier.hh"

#include "cuda_utils/Shared_Device_Ptr.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Domain_Transporter
 * \brief Transport a particle on a computational domain.
 *
 * This class does no communication; it takes a particle and transports it
 * until it leaves the domain.
 */
/*!
 * \example cuda_mc/test/tstDomain_Transporter.cc
 *
 * Test of Domain_Transporter.
 */
//===========================================================================//

template <class Geometry>
class Domain_Transporter
{
  public:
    //@{
    //! Useful typedefs.
    typedef Geometry					Geometry_t;
    typedef Physics<Geometry_t>				Physics_t;
    typedef typename Geometry_t::Space_Vector		Space_Vector;
    typedef typename Geometry_t::Geo_State_t		Geo_State_t;
    typedef typename Physics_t::Particle_Vector_t       Particle_Vector_t;
    typedef typename Physics_t::Bank_t			Bank_t;
    typedef typename Physics_t::Fission_Site_Container	Fission_Site_Container;
    typedef Variance_Reduction<Geometry_t>		Variance_Reduction_t;
    typedef Tallier<Geometry_t>				Tallier_t;
    //@}

    //@{
    //! Smart pointers.
    typedef cuda::Shared_Device_Ptr<Geometry_t>		SDP_Geometry;
    typedef cuda::Shared_Device_Ptr<Physics_t>		SDP_Physics;
    typedef cuda::Shared_Device_Ptr<Particle_Vector_t>	SDP_Particle_Vector;
    typedef cuda::Shared_Device_Ptr<Bank_t>		SDP_Bank;
    typedef std::shared_ptr<Variance_Reduction_t>	SP_Variance_Reduction;
    typedef std::shared_ptr<Tallier_t>			SP_Tallier;
    typedef std::shared_ptr<Fission_Site_Container>	SP_Fission_Sites;
    //@}

  private:
    // >>> DATA

    // Problem geometry implementation.
    SDP_Geometry d_geometry;

    // Problem physics implementation.
    SDP_Physics d_physics;

    // Variance reduction.
    SP_Variance_Reduction d_var_reduction;

    // Regular tallies.
    SP_Tallier d_tallier;

    // Fission sites.
    SP_Fission_Sites d_fission_sites;

  public:
    // Constructor.
    Domain_Transporter();

    // Set the geometry and physics classes.
    void set(const SDP_Geometry& geometry, const SDP_Physics& physics);

    // Set the variance reduction.
    void set(const SP_Variance_Reduction& reduction);

    // Set regular tallies.
    void set(const SP_Tallier& tallies);

    // Set fission site sampling.
    void set(const SP_Fission_Sites& fission_sites, double keff);

    // Transport a particle one step through the domain to either a collision
    // or a boundary and do the pathlength tallies.
    void transport_step(SDP_Particle_Vector_t& particles, SDP_Bank_t& bank);

    // Process particles that have hit a boundary.
    void process_boundary(SDP_Particle_Vector_t& particles, SDP_Bank_t& bank);

    // Process particles that have hit a collision.
    void process_collision(SDP_Particle_Vector_t& particles, SDP_Bank_t& bank);

    //! Return the number of sampled fission sites.
    int num_sampled_fission_sites() const { return d_num_fission_sites; }

  private:
    // >>> IMPLEMENTATION

    // Flag indicating that fission sites should be sampled.
    bool d_sample_fission_sites;

    // Number of fission sites sampled.
    int d_num_fission_sites;

    // Current keff iterate.
    double d_keff;

    // Process collisions and boundaries.
    void process_boundary(SDP_Particle_Vector_t &particle, Bank_t &bank);
    void process_collision(SDP_Particle_Vector_t &particle, Bank_t &bank);
};

} // end namespace cuda_profugus

#endif // cuda_mc_Domain_Transporter_hh

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.hh
//---------------------------------------------------------------------------//
