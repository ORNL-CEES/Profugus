//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Source.cuh
 * \author Steven Hamilton
 * \date   Mon May 05 14:22:41 2014
 * \brief  Source class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Source_cuh
#define cuda_mc_Source_cuh

#include <iostream>
#include <memory>

#include "utils/Definitions.hh"
#include "cuda_utils/Definitions.hh"
#include "cuda_utils/Shared_Device_Ptr.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Source
 * \brief Base class definition for Monte Carlo sources.
 */
//===========================================================================//

template <class Geometry>
class Source
{
  public:
    //@{
    //! Useful typedefs.
    typedef Geometry                            Geometry_t;
    typedef typename Geometry_t::Geo_State_t    Geo_State_t;
    typedef cuda_profugus::Space_Vector         Space_Vector;
    typedef def::size_type                      size_type;
    //@}

    //@{
    //! Smart pointers
    typedef cuda::Shared_Device_Ptr<Geometry_t> SDP_Geometry;
    //@}

  protected:
    // >>> DATA

    // Geometry and physics.
    Geometry_t  *b_geometry;

    // Requested particles per cycle.
    size_type d_np_requested;

    // Number of particles: total, domain
    size_type d_np_total, d_np_domain;
    
    // Number of particles run/left on this domain
    size_type d_np_run;
    size_type d_np_left;

    // Size of batch
    size_type d_batch_size;

  public:
    // Constructor.
    Source(SDP_Geometry geometry);

    // Virtual destructor for polymorphism.
    virtual ~Source(){}

    // Derived classes should implement the following functions,
    // but because this is an on-device class there can be NO virtual
    // functions.  Toggling between derived classes must be handled
    // through templates.

    //! Get a particle from the source on specified thread.
    //__device__ Particle_t get_particle(
    //    std::size_t idx, RNG_State_t &rng) const;

    //! Number of particles to transport in the source on the current domain.
    size_type num_to_transport() const {return d_np_domain;}

    //! Total number of particles to transport in the entire problem/cycle.
    size_type total_num_to_transport() const {return d_np_total;}

    //! Number of particles remaining to transport on current domain
    size_type num_left() const {return d_np_left;}

    //! Number of particles to transport in current batch
    size_type num_batch() const
    {
        REQUIRE(d_batch_size > 0);
        return std::min(d_np_left,d_batch_size);
    }

    //! Prepare for new batch of source particles
    virtual void begin_batch(size_type Np){};

    //! Set batch size
    virtual void set_batch_size(size_type batch_size)
    {
        REQUIRE(batch_size>0);
        d_batch_size = std::min(batch_size,d_np_domain);
    }

    //! Update number of histories left after completing batch
    void end_batch(size_type Np)
    {
        REQUIRE(Np <= d_np_left);
        d_np_left -= Np;
        d_np_run  += Np;
    }

    // >>> INHERITED INTERFACE

    //! Get the geometry.
    const Geometry_t& geometry() const { return *b_geometry; }
};

} // end namespace cuda_mc

#endif // cuda_mc_Source_cuh

//---------------------------------------------------------------------------//
//                 end of Source.cuh
//---------------------------------------------------------------------------//
