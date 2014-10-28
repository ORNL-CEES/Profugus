//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/mc/Sampler.hh
 * \author Thomas M. Evans
 * \date   Friday May 2 10:26:10 2014
 * \brief  Sampling functions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef core_mc_Sampler_hh
#define core_mc_Sampler_hh

namespace profugus
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
int sample_discrete_CDF(int nb, const T *c, const T ran);

// Sample a discrete CDF using iterators and a logarithmic search.
template<class InputIterator, class FloatType>
inline InputIterator sample_dcdf(
        FloatType      xi,
        InputIterator  iter,
        InputIterator  last);

// Sample a discrete CDF using iterators and a linear search.
template<class InputIterator, class FloatType>
InputIterator sample_small_dcdf(
        FloatType      xi,
        InputIterator  first,
        InputIterator  last);

// Sample a discrete CDF using a linear search starting at the back
template<class InputIterator, class FloatType>
InputIterator sample_smallrev_dcdf(
        FloatType      xi,
        InputIterator  first,
        InputIterator  last);

// Sample an energy in eV from a Watt fission spectrum.
template<class RNG>
double sample_watt(RNG &rng, double a = 0.965, double b = 2.29);

// Sample normalized linear distribution on [0,1)
template<class T>
inline T sample_linear(T xi);

// Sample a linear distribution on [0,1) with arbitrary y values
template<class T, class RNG>
inline T sample_linear(RNG& rng, const T left, const T right);

// Sample a linear distribution on [0,1) with arbitrary y values
template<class T>
T sample_linear(const T xi_lr, const T xi_linear, const T left, const T right);

//---------------------------------------------------------------------------//
} // end namespace sampler

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "Sampler.i.hh"

//---------------------------------------------------------------------------//
#endif // core_mc_Sampler_hh

//---------------------------------------------------------------------------//
//              end of mc/Sampler.hh
//---------------------------------------------------------------------------//
