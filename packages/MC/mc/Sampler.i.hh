//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Sampler.i.hh
 * \author Seth R Johnson
 * \date   Friday May 2 10:26:27 2014
 * \brief  Member definitions of class Sampler.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Sampler_i_hh
#define mc_Sampler_i_hh

#include <algorithm>
#include <cmath>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"

namespace profugus
{

namespace sampler
{

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a discrete CDF using iterators and a logarithmic search.
 *
 * This assumes a uniform input in the range
 *
 * \param iter Beginning of DCDF
 * \param last End of DCDF
 * \param xi A random number uniformly generated in [0,1)
 */
template<class InputIterator, class FloatType>
inline InputIterator sample_dcdf(
        FloatType      xi,
        InputIterator  iter,
        InputIterator  last)
{
    Require(std::distance(iter, last) > 0);
    Require(xi >= 0 && xi < 1);

    // Do a binary search on the CDF
    iter = std::lower_bound(iter, last, xi);

    // Return the value
    Ensure(iter != last);
    return iter;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a discrete CDF using iterators and a linear search.
 *
 * \param iter Beginning of DCDF
 * \param last End of DCDF
 * \param xi A random number uniformly generated in [0,1)
 */
template<class InputIterator, class FloatType>
inline InputIterator sample_small_dcdf(
        FloatType      xi,
        InputIterator  iter,
        InputIterator  last)
{
    Require(std::distance(iter, last) > 0);
    Require(xi >= 0 && xi < 1);

    while (iter != last)
    {
        if (*iter >= xi)
            break;
        ++iter;
    }
    return iter;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a discrete CDF using a reversed linear search
 *
 * This will reduce the sampling time when elements are more likely to be
 * toward the back of the distribution.
 *
 * \param iter Beginning of DCDF
 * \param last End of DCDF
 * \param xi A random number uniformly generated in [0,1)
 */
template<class InputIterator, class FloatType>
inline InputIterator sample_smallrev_dcdf(
        FloatType      xi,
        InputIterator  first,
        InputIterator  iter)
{
    Require(std::distance(first, iter) > 0);
    Require(xi >= 0 && xi < 1);

    while (iter != first)
    {
        --iter;
        if (*iter < xi)
        {
            ++iter;
            break;
        }
    }
    return iter;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a normalized linear distribution on [0,1)
 */
template<class T>
inline T sample_linear(T xi)
{
    Require(0 <= xi && xi < 1);

    return std::sqrt(xi);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a linear distribution on [0,1) with arbitrary y values
 *
 * This pulls two random numbers from the generator and passes them to the
 * underlying sampler function.
 */
template<class T, class RNG>
inline T sample_linear(RNG& rng, const T left, const T right)
{
    return sample_linear(
            rng.template uniform<T>(), rng.template uniform<T>(),
            left, right);
}

//---------------------------------------------------------------------------//

} // end namespace profugus::sampler

} // end namespace profugus

#endif // mc_Sampler_i_hh

//---------------------------------------------------------------------------//
//                 end of Sampler.i.hh
//---------------------------------------------------------------------------//
