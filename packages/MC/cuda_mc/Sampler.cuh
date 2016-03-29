//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Sampler.cuh
 * \author Steven Hamilton
 * \date   Friday May 2 10:26:10 2014
 * \brief  Sampling functions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Sampler_cuh
#define cuda_mc_Sampler_cuh

#include "Definitions.cuh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
/*!
 * \example mc/test/tstSampler.cc
 *
 * Test of sampler namespace functions.
 */
//---------------------------------------------------------------------------//

//! Sampling functions namespace
namespace sampler
{

//---------------------------------------------------------------------------//
// Sample a discrete CDF.
template<class T>
__device__ inline int sample_discrete_CDF(int nb, const T *c, const T ran);

// Sample a discrete CDF using iterators and a logarithmic search.
template<class FloatType>
__device__ inline FloatType * sample_dcdf(
        FloatType   xi,
        FloatType  *iter,
        FloatType  *last);

// Sample a discrete CDF using iterators and a linear search.
template<class FloatType>
__device__ inline FloatType * sample_small_dcdf(
        FloatType   xi,
        FloatType  *first,
        FloatType  *last);

// Sample a discrete CDF using a linear search starting at the back
template<class FloatType>
__device__ inline FloatType * sample_smallrev_dcdf(
        FloatType   xi,
        FloatType  *first,
        FloatType  *last);

// Sample an energy in eV from a Watt fission spectrum.
template<class RNG_State>
__device__ inline double sample_watt(RNG_State *rng,
                                     double     a = 0.965,
                                     double     b = 2.29);

// Sample normalized linear distribution on [0,1)
template<class T>
__device__ inline T sample_linear(T xi);

// Sample a linear distribution on [0,1) with arbitrary y values
template<class RNG_State>
__device__ inline float sample_linear(RNG_State *rng, const float left,
                                      const float right);
template<class RNG_State>
__device__ inline double sample_linear(RNG_State *rng, const double left,
                                       const double right);

// Sample a linear distribution on [0,1) with arbitrary y values
template<class T>
__device__ inline T sample_linear(const T xi_lr, const T xi_linear,
                                  const T left,  const T right);

// Sample an isotropic angular distribution
__device__ inline void sample_isotropic(cuda::Space_Vector &omega,
                                        RNG_State_t        *rng);

//---------------------------------------------------------------------------//
} // end namespace sampler

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "Sampler.i.cuh"

//---------------------------------------------------------------------------//
#endif // cuda_mc_Sampler_cuh

//---------------------------------------------------------------------------//
//              end of cuda_mc/Sampler.cuh
//---------------------------------------------------------------------------//
