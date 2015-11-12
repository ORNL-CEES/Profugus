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
#include "Uniform_Source.hh"
#include "config.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>

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
    typedef Uniform_Source                             Source_t;
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
    typedef std::shared_ptr<Source_t>               SP_Source;
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

    // Particle source.
    SP_Source d_source;

  public:
    // Constructor.
    Domain_Transporter();

    // Set the geometry and physics classes.
    void set(SP_Geometry geometry, SP_Physics physics);

    // Set the variance reduction.
    void set(SP_Variance_Reduction reduction);

    // Set regular tallies.
    void set(SP_Tallier tallies);

    // Set the source.
    void set(SP_Source source);

    // Transport a particle through the domain.
    void transport( std::vector<Particle_t>& particles,
		    std::vector<std::pair<std::size_t,events::Event> >& events, 
		    std::vector<Bank_t>& banks ) const;

  private:
    // >>> IMPLEMENTATION

    // Process collisions and boundaries.
    void launch_cuda() const;
    double distance_to_collision( Particle_t& particle,
				  const double& xs_total ) const;
    void get_next_event( Particle_t& particle, events::Event& event ) const;
    void process_boundary(
	Particle_t &particle, events::Event& event, Bank_t &bank) const;
    void process_collision(
	Particle_t &particle, events::Event& event, Bank_t &bank) const;

  private:

#ifdef USE_TRILINOS_TIMING
    // Get next event timer.
    Teuchos::RCP<Teuchos::Time> d_next_event_timer;

    // Get event sort timer.
    Teuchos::RCP<Teuchos::Time> d_event_sort_timer;

    // Process collision timer.
    Teuchos::RCP<Teuchos::Time> d_process_collision_timer;

    // Process boundary timer.
    Teuchos::RCP<Teuchos::Time> d_process_boundary_timer;
    
    // New mfp timer.
    Teuchos::RCP<Teuchos::Time> d_new_mfp_timer;

    // Tally timer.
    Teuchos::RCP<Teuchos::Time> d_tally_timer;    

    // Alive sort timer.
    Teuchos::RCP<Teuchos::Time> d_alive_sort_timer;    
#endif
    
};

} // end namespace profugus

#endif // mc_Domain_Transporter_hh

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.hh
//---------------------------------------------------------------------------//
