//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Source.hh
 * \author Thomas M. Evans
 * \date   Mon May 05 14:22:41 2014
 * \brief  Source class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Source_hh
#define mc_Source_hh

#include <memory>
#include <cmath>

#include "utils/Definitions.hh"
#include "geometry/Geometry.hh"
#include "Physics.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Source
 * \brief Base class definition for Monte Carlo sources.
 */
//===========================================================================//

class Source
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

  protected:
    // >>> DATA

    // Geometry and physics.
    SP_Geometry b_geometry;
    SP_Physics  b_physics;

    // Sample isotropic angle.
    void sample_angle(Space_Vector &omega, RNG& rng)
    {
        using def::X; using def::Y; using def::Z;

        omega[Z]        = 1.0 - 2.0 * rng.ran();
        double phi      = constants::two_pi * rng.ran();
        double sintheta = std::sqrt(1.0 - omega[Z] * omega[Z]);

        omega[X] = sintheta * std::cos(phi);
        omega[Y] = sintheta * std::sin(phi);
    }

    // Node ids.
    const int b_node, b_nodes;

  public:
    // Constructor.
    Source(SP_Geometry    geometry,
           SP_Physics     physics );

    //! Virtual destructor.
    virtual ~Source() = 0;

    // >>> VIRTUAL PUBLIC INTERFACE

    //! Get a particle from the source.
    virtual SP_Particle get_particle() = 0;

    //! Whether the source has finished emitting all its particles.
    virtual bool empty() const = 0;

    //! Number of particles to transport in the source on the current domain.
    virtual size_type num_to_transport() const = 0;

    //! Total number of particles to transport in the entire problem/cycle.
    virtual size_type total_num_to_transport() const = 0;

    // >>> INHERITED INTERFACE

    //! Get the geometry.
    const Geometry_t& geometry() const { return *b_geometry; }

    //! Get the physics.
    const Physics_t& physics() const { return *b_physics; }
};

} // end namespace profugus

#endif // mc_Source_hh

//---------------------------------------------------------------------------//
//                 end of Source.hh
//---------------------------------------------------------------------------//
