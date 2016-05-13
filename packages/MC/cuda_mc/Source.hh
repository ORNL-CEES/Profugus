//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Source.hh
 * \author Stuart Slattery
 * \date   Mon May 05 14:22:41 2014
 * \brief  Source class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Source_hh
#define cuda_mc_Source_hh

#include "cuda_utils/Shared_Device_Ptr.hh"

#include "Definitions.hh"
#include "Particle_Vector.hh"

namespace cuda_profugus
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

    //! Virtual destructor.
    virtual ~Source() { /* ... */ };

    // >>> VIRTUAL PUBLIC INTERFACE

    //! Get particles from the source.
    virtual void get_particles( 
	cuda::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles ) = 0;

    //! Whether the source has finished emitting all its particles.
    virtual bool empty() const = 0;

    //! Number of particles to transport in the source on the current domain.
    virtual std::size_t num_to_transport() const = 0;

    //! Total number of particles to transport in the entire problem/cycle.
    virtual std::size_t total_num_to_transport() const = 0;

    //! Total number of requested particles.
    virtual std::size_t Np() const = 0;

    //! Number transported so far on this domain.
    virtual std::size_t num_run() const = 0;

    //! Number left to transport on this domain.
    virtual std::size_t num_left() const = 0;
};

} // end namespace cuda_profugus

#endif // cuda_mc_Source_hh

//---------------------------------------------------------------------------//
//                 end of Source.hh
//---------------------------------------------------------------------------//
