//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Particle.hh
 * \author Thomas M. Evans
 * \date   Fri Apr 25 11:26:16 2014
 * \brief  Particle class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Particle_hh
#define cuda_mc_Particle_hh

#include "utils/Definitions.hh"
#include "cuda_utils/CudaMacros.hh"
#include "rng/RNG.hh"
#include "geometry/RTK_State.hh"
#include "Definitions.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Particle
 * \brief Particle class for MC transport.
 */
/*!
 * \example mc/test/tstParticle.cc
 *
 * Test of Particle.
 */
//===========================================================================//
template <class Geometry>
class Particle
{
  public:
    //@{
    //! Typedefs.
    typedef def::Space_Vector               Space_Vector;
    typedef events::Event                   Event_Type;
    typedef typename Geometry::Geo_State_t  Geo_State_t;
    //@}

  private:
    // >>> DATA

    // Material id in current region.
    int d_matid;

    // Particle group index.
    int d_group;

    // Particle weight.
    double d_wt;

    // Alive/dead status.
    bool d_alive;

    // Latest particle event.
    Event_Type d_event;

    // Particle geometric state.
    Geo_State_t d_geo_state;

    // Particle batch.
    int d_batch;

  public:
    // >>> PARTICLE FUNCTIONS

    //! Set a new weight.
    PROFUGUS_HOST_DEVICE_FUNCTION
    void set_wt(double wt) { d_wt = wt; }

    //! Set particle group.
    PROFUGUS_HOST_DEVICE_FUNCTION
    void set_group(int g) { d_group = g; }

    //! Multiply weight.
    PROFUGUS_HOST_DEVICE_FUNCTION
    void multiply_wt(double wt) { d_wt *= wt; }

    //! Set the particle event flag.
    PROFUGUS_HOST_DEVICE_FUNCTION
    void set_event(Event_Type event) { d_event = event; }

    //! Set the material id of the region occupied by the particle.
    PROFUGUS_HOST_DEVICE_FUNCTION
    void set_matid(int matid) { d_matid = matid; }

    //! Kill the particle.
    PROFUGUS_HOST_DEVICE_FUNCTION
    void kill() { d_alive = false; }

    //! Set particle status to alive.
    PROFUGUS_HOST_DEVICE_FUNCTION
    void live() { d_alive = true; }

    //! Set the particle batch.
    PROFUGUS_HOST_DEVICE_FUNCTION
    void set_batch(int b) { d_batch = b; }

    //@{
    //! Get a handle to the geometric state of the particle.
    PROFUGUS_HOST_DEVICE_FUNCTION
    Geo_State_t& geo_state() { return d_geo_state; }

    PROFUGUS_HOST_DEVICE_FUNCTION
    const Geo_State_t& geo_state() const { return d_geo_state; }
    //@}

    //@{
    //! Access particle data.
    PROFUGUS_HOST_DEVICE_FUNCTION
    bool alive() const { return d_alive; }

    PROFUGUS_HOST_DEVICE_FUNCTION
    double wt() const { return d_wt; }

    PROFUGUS_HOST_DEVICE_FUNCTION
    Event_Type event() const { return d_event; }

    PROFUGUS_HOST_DEVICE_FUNCTION
    int matid() const { return d_matid; }

    PROFUGUS_HOST_DEVICE_FUNCTION
    int group() const { return d_group; }

    PROFUGUS_HOST_DEVICE_FUNCTION
    int batch() const{ return d_batch; }
    //@}
};

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Particle_hh

//---------------------------------------------------------------------------//
//                 end of Particle.hh
//---------------------------------------------------------------------------//
