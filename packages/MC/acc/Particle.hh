//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   acc/Particle.hh
 * \author Thomas M. Evans
 * \date   Tue Oct 28 12:40:21 2014
 * \brief  OpenACC Particle implementation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef acc_Particle_hh
#define acc_Particle_hh

#include <vector>

#include "core/geometry/Mesh_State.hh"
#include "core/geometry/Mesh_Geometry.hh"
#include "core/mc/Definitions.hh"
#include "core/mc/Source.hh"

namespace acc
{

//===========================================================================//
/*!
 * \class Particle
 * \brief ACC Particle type.
 */
//===========================================================================//

struct Particle
{
    typedef profugus::Mesh_State    Geo_State_t;
    typedef profugus::events::Event Event_Type;

    // >>> DATA

    // Material id in current region.
    int matid;

    // Particle group index.
    int group;

    // Particle weight.
    double wt;

    // Alive/dead status.
    bool alive;

    // Latest particle event.
    Event_Type event;

    // Particle geometric state.
    Geo_State_t geo_state;
};

//! ACC particle containers.
typedef std::vector<Particle> Vec_Particles;

//! Random numbers.
extern std::vector<double> rnd_numbers;

//! Set the ideal vector size of particles.
void set_size(Vec_Particles &particles, int num_steps);

//! Load particles from a source.
void load_source(profugus::Mesh_Geometry &geometry, profugus::Source &source,
                 Vec_Particles &particles);

} // end namespace acc

#endif // acc_Particle_hh

//---------------------------------------------------------------------------//
//                 end of Particle.hh
//---------------------------------------------------------------------------//
