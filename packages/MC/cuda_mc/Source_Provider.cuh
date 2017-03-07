//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/Source_Provider.cuh
 * \author Steven Hamilton
 * \date   Wed Apr 13 08:38:06 2016
 * \brief  Source_Provider class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_mc_Source_Provider_cuh
#define MC_cuda_mc_Source_Provider_cuh

#include <thrust/device_vector.h>
#include "Definitions.hh"
#include "Source.cuh"
#include "Particle_Vector.cuh"
#include "RNG_Control.cuh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Source_Provider
 * \brief Compute vector of particles from a given source.
 *
 * Given a source object (which implements the on-device method "get_particle"
 * this class computes a thrust::device_vector of Particles representing
 * the initial states of the source particles.
 *
 * This class provides an opportunity to manipulate the initial order of
 * particles (e.g. sorted by cell or material id) to introduce some locality
 * in data accesses.  It can also control the batching of histories (for
 * limiting memory usage and/or computing batch statistics).  
 *
 * \sa Source_Provider.t.cuh for detailed descriptions.
 */
/*!
 * \example cuda_mc/test/tstSource_Provider.cc
 *
 * Test of Source_Provider.
 */
//===========================================================================//

template <class Geometry>
class Source_Provider
{
  public:

    typedef Particle_Vector<Geometry>               Particle_Vector_t;
    typedef Particle_Vector_DMM<Geometry>           Particle_Vector_DMM_t;
    typedef std::shared_ptr<Particle_Vector_DMM_t>  SP_Particle_Vector_DMM;
    typedef int                                     index;
    typedef thrust::device_vector<index>            Index_Vector;
    typedef std::shared_ptr<RNG_Control>            SP_RNG_Control;
    typedef std::shared_ptr<Source<Geometry>>       SP_Source;
    typedef def::size_type                          size_type;

    //! Constructor
    Source_Provider(){}

    //! Get vector of particles
    void get_particles( SP_Source               source,
                        SP_RNG_Control          rng_control,
                        SP_Particle_Vector_DMM  particles,
                        const Index_Vector     &indices ) const;

  private:

    // Implementation of get_particles for particular source type
    template <class Src_Type>
    void get_particles_impl( Src_Type                  &source,
                             SP_RNG_Control             rng_control,
                             SP_Particle_Vector_DMM     particles,
                             const Index_Vector        &indices,
                             size_type                  Np) const;
                        
};

//---------------------------------------------------------------------------//
} // end namespace cuda_mc

#endif // MC_cuda_mc_Source_Provider_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Source_Provider.cuh
//---------------------------------------------------------------------------//
