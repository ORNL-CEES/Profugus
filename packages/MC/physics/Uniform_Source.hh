//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Uniform_Source.hh
 * \author Thomas M. Evans
 * \date   Tue May 06 16:43:26 2014
 * \brief  Uniform_Source class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Uniform_Source_hh
#define mc_Uniform_Source_hh

#include <hpx/include/components.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/util.hpp>

#include <vector>
#include <memory>
#include <cmath>

#include "utils/Definitions.hh"
#include "geometry/Geometry.hh"
#include "Physics.hh"
#include "Shape.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Uniform_Source
 * \brief Defines a uniformly sampled source for fixed-source problems.
 *
 * Currently, the only implemented source shape is a box (see Box_Shape).
 * Also, all angular sampling is isotropic.
 *
 * \section uniform_source_db Standard DB Entries for Uniform_Source
 *
 * The following standard data entries control the Uniform_Source:
 *
 * \arg \c Np (int) number of particles to use in each cycle (default:
 * 1000)
 *
 * \arg \c spectral_shape (Array<double>) source spectral (energy) shape by
 * group (default: flat)
 */
/*!
 * \example mc/test/tstUniform_Source.cc
 *
 * Test of Uniform_Source.
 */
//===========================================================================//
class Uniform_Source : public hpx::components::locking_hook<
    hpx::components::managed_component_base<Uniform_Source> >
{
  public:
    //@{
    //! Useful typedefs.
    typedef Core                     Geometry_t;
    typedef Geometry_t::Geo_State_t  Geo_State_t;
    typedef Physics                  Physics_t;
    typedef Physics_t::Particle_t    Particle_t;
    typedef Geometry_t::Space_Vector Space_Vector;
    typedef def::size_type           size_type;
    //@}

    //@{
    //! Smart pointers
    typedef std::shared_ptr<Geometry_t>    SP_Geometry;
    typedef std::shared_ptr<Physics_t>     SP_Physics;
    typedef std::shared_ptr<Particle_t>    SP_Particle;
    //@}

    //@{
    //! Typedefs.
    typedef Physics_t::RCP_Std_DB  RCP_Std_DB;
    typedef def::Vec_Dbl           Vec_Dbl;
    typedef std::shared_ptr<Shape> SP_Shape;
    //@}

  private:
    // >>> DATA

    // Geometry and physics.
    SP_Geometry d_geometry;
    SP_Physics  d_physics;

    // Source geometric shape.
    SP_Shape d_geo_shape;

    // Energy shape CDF.
    Vec_Dbl d_erg_cdf;

    // Node ids.
    const int d_node, d_nodes;

  public:
    // Default constructor for HPX.
    Uniform_Source();

    // Constructor.
    Uniform_Source(RCP_Std_DB db, SP_Geometry geometry, SP_Physics physics, SP_Shape geometric_shape);

    // >>> DERIVED PUBLIC INTERFACE

    // Get a particle from the source.
    SP_Particle get_particle( const int lid );

    //! Number of particles to transport in the source on the current domain.
    size_type num_to_transport() const { return d_np_domain; }

    //! Total number of particles to transport in the entire problem/cycle.
    size_type total_num_to_transport() const { return d_np_total; }

    //! Get the geometry.
    const Geometry_t& geometry() const { return *d_geometry; }

    //! Get the physics.
    const Physics_t& physics() const { return *d_physics; }

    // >>> CLASS ACCESSORS

    //! Total number of requested particles.
    size_type Np() const { return d_np_requested; }

    // >>> HPX DEFINITIONS

    // Register the get particle function as an action.
    HPX_DEFINE_COMPONENT_ACTION( Uniform_Source, get_particle );

  private:

    // Build the domain replicated source.
    void build_DR();

    // Sample isotropic angle.
    static Space_Vector sample_angle(const double r1, const double r2)
    {
        using def::X; using def::Y; using def::Z;

	Space_Vector omega;

        omega[Z]        = 1.0 - 2.0 * r1;
        double phi      = constants::two_pi * r2;
        double sintheta = std::sqrt(1.0 - omega[Z] * omega[Z]);

        omega[X] = sintheta * std::cos(phi);
        omega[Y] = sintheta * std::sin(phi);

	return omega;
    }

  private:
    // >>> IMPLEMENTATION

    // Requested particles.
    size_type d_np_requested;

    // Number of particles: total, domain
    size_type d_np_total;
    size_type d_np_domain;

    // Particle weight.
    const double d_wt;
};

} // end namespace profugus

//---------------------------------------------------------------------------//
// HPX GLOBAL NAMESPACE DECLARATION
//---------------------------------------------------------------------------//
// Register the get_particle function.
HPX_REGISTER_ACTION_DECLARATION( profugus::Uniform_Source::get_particle_action,
				 uniform_source_get_particle_action );

//---------------------------------------------------------------------------//
// HPX CLIENT CLASS
//---------------------------------------------------------------------------//
namespace profugus
{
class Uniform_Source_Client : public hpx::components::client_base<
    Uniform_Source_Client,Uniform_Source>
{
  public:

    using base_type =
	hpx::components::client_base<Uniform_Source_Client,Uniform_Source>;

    // Default constructor.
    Uniform_Source_Client()
    { /* ... */ }

    // Create a client for an existing server with a given GID.
    Uniform_Source_Client( hpx::future<hpx::naming::id_type>&& gid )
	: base_type( std::move(gid) )
    { /* ... */ }

    // Get a particle as an asynchronous task.
    hpx::future<Uniform_Source::SP_Particle> get_particle_async( const int lid )
    {
	HPX_ASSERT( this->get_id() );
	return hpx::async<Uniform_Source::get_particle_action>( 
	    this->get_id(), lid );
    }
};

} // and namespace profugus

//---------------------------------------------------------------------------//

#endif // mc_Uniform_Source_hh

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.hh
//---------------------------------------------------------------------------//
