//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Domain_Transporter.hh
 * \author Thomas M. Evans
 * \date   Mon May 12 12:02:13 2014
 * \brief  Domain_Transporter class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Domain_Transporter_hh
#define mc_Domain_Transporter_hh

#include <memory>

#include "geometry/Geometry.hh"
#include "Physics.hh"
#include "Step_Selector.hh"
#include "Variance_Reduction.hh"
#include "Tallier.hh"

namespace profugus
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
 * \example mc/test/tstDomain_Transporter.cc
 *
 * Test of Domain_Transporter.
 */
//===========================================================================//

class Domain_Transporter
{
  public:
    //@{
    //! Useful typedefs.
    typedef Core                                       Geometry_t;
    typedef Physics                                    Physics_t;
    typedef typename Geometry_t::Space_Vector          Space_Vector;
    typedef typename Geometry_t::Geo_State_t           Geo_State_t;
    typedef typename Physics_t::Particle_t             Particle_t;
    typedef typename Physics_t::Bank_t                 Bank_t;
    typedef Variance_Reduction                         Variance_Reduction_t;
    typedef Tallier                                    Tallier_t;
    //@}

    //@{
    //! Smart pointers.
    typedef std::shared_ptr<Geometry_t>             SP_Geometry;
    typedef std::shared_ptr<Physics_t>              SP_Physics;
    typedef std::shared_ptr<Variance_Reduction_t>   SP_Variance_Reduction;
    typedef std::shared_ptr<Tallier_t>              SP_Tallier;
    //@}

  private:
    // >>> DATA

    // Problem geometry implementation.
    SP_Geometry d_geometry;

    // Problem physics implementation.
    SP_Physics d_physics;

    // Variance reduction.
    SP_Variance_Reduction d_var_reduction;

    // Regular tallies.
    SP_Tallier d_tallier;

  public:
    // Constructor.
    Domain_Transporter();

    // Set the geometry and physics classes.
    void set(SP_Geometry geometry, SP_Physics physics);

    // Set the variance reduction.
    void set(SP_Variance_Reduction reduction);

    // Set regular tallies.
    void set(SP_Tallier tallies);

    // Transport a particle through the domain.
    void transport(Particle_t &particle, 
		   events::Event& event, 
		   Bank_t &bank) const;

  private:
    // >>> IMPLEMENTATION

    // Process collisions and boundaries.
    void get_next_event( Particle_t& particle, events::Event& event, double& dist_mfp ) const;
    void process_boundary(Particle_t &particle, events::Event& event, Bank_t &bank) const;
    void process_collision(
	const double step, Particle_t &particle, events::Event& event, Bank_t &bank) const;
};

} // end namespace profugus

#endif // mc_Domain_Transporter_hh

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.hh
//---------------------------------------------------------------------------//
