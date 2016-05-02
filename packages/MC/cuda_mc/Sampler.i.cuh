//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Sampler.i.cuh
 * \author Steven Hamilton
 * \date   Friday May 2 10:26:27 2014
 * \brief  Member definitions of class Sampler.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Sampler_i_cuh
#define cuda_mc_Sampler_i_cuh

#include "cuda_utils/Utility_Functions.hh"
#include "cuda_utils/Constants.hh"

namespace cuda_mc
{

namespace sampler
{

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a discrete CDF.
 *
 * Given a CDF defined over the range \f$[0,N_b)\f$, where \f$N_b\f$ is the * number of bins, sample and return the bin of the sampled value. The domain
 * of the CDF is \f$[0,1]\f$.  The sampling works as follows,
 * \verbatim
   Bin    0     1     2      3      4
       |     |     |     |       |     |
       |-----|-----|-----|-------|-----|
       |     |     |     |       |     |
      0.0    c[0]  c[1]  c[2]    c[3]  c[4] = 1.0
 * \endverbatim
 * Here, \f$N_b=5\f$.
 *
 * \param nb       number of bins
 * \param c        pointer of type T (float or double) pointing to CDF
 * \param ran      random number between \f$[0.0,1.0]\f$
 * \return sampled bin index between \f$[0,N_b)\f$
 */

template<class T>
__device__  int sample_discrete_CDF(int      nb,
                                    const T *c,
                                    const T  ran)
{
    REQUIRE(nb > 0);
    REQUIRE(cuda::utility::soft_equiv(static_cast<double>(c[nb - 1]), 1.0,
                                      1.0e-6));
    REQUIRE(ran >= 0 && ran <= 1);

    // do a binary search on the CDF
    const T *ptr = cuda::utility::lower_bound(c, c + nb, ran);

    // return the value
    ENSURE(ptr - c >= 0 && ptr - c < nb);
    return ptr - c;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample an energy in eV from a Watt fission spectrum.
 *
 * Given a reference to a random number generator and constants \f$a\f$ and
 * \f$b\f$, generate a sample from the Watt fission spectrum:
 * \f[
 *     f(E) = c * exp(-E/a) * sinh(sqrt(b*E))
 * \f]
 * where \f$c\f$ is a normalization constant. The default values of \f$a\f$ and
 * \f$b\f$ here are the same as the defaults used in MCNP5. The random number
 * generator must provide a \c ran() function returning a double in the range
 * \f$(0,1]\f$.
 *
 * \param rng reference of type T to a random number generator
 * \param a   first constant (in MeV)
 * \param b   second constant (in 1/MeV)
 * \return sampled energy in \f$(0,\infty)\f$ in eV
 */
template <class RNG_State>
__device__  double sample_watt(RNG_State   *rng,
                               double       a,
                               double       b)
{
    REQUIRE(a > 0);
    REQUIRE(b > 0);

    // the algorithm implemented here is the same as in MCNP5

    const double t1 = 1 + 0.125 * a * b;
    const double c1 = t1 + sqrt(t1 * t1 - 1) - 1;
    const double c2 = a * (c1 + 1);
    const double c3 = b * c2;

    while (1)
    {
        double r1 = -log(curand_uniform_double(rng));
        double r2 = -log(curand_uniform_double(rng));
        double c4 = r2 - c1 * (1 + r1);
        if (c4 * c4 <= c3 * r1)
            return 1.0e6 * c2 * r1;  // in eV
    }
}


//---------------------------------------------------------------------------//
/*!
 * \brief Sample a linear distribution on [0,1) with arbitrary y values.
 *
 * This sampling function requires two random numbers. Because of how Scale CE
 * physics uses random numbers in its interpolation schemes, we provide this
 * function explicitly as a function of those two numbers.
 *
 * Source:
 * Forrest B. Brown, Fundamentals of Monte Carlo particle transport,
 * Tech. Report LA-UR-04-8817, Los Alamos National Laboratory, 2004.
 *
 * \param xi_lr Random number for whether to sample the left or right triangle
 * \param xi_linear Random number for sampling inside the triangle
 * \param left Left linear PDF value
 * \param right Right linear PDF value
 *
 * \return sampled value on [0,1)
 */
template<class T>
__device__  T sample_linear(
        const T xi_lr,
        const T xi_linear,
        const T left,
        const T right)
{
    REQUIRE(0 <= xi_lr && xi_lr <= 1);
    REQUIRE(0 <= xi_linear && xi_linear <= 1);

    // use the method that decomposes the line into basically two linear
    // discontinous basis functions
    // total area: (1 - 0) * (rightval - leftval)/2 + leftval * (1 - 0)
    //               = rightval/2 + leftval/2 (same as two triangles)
    // probability of being in the left: leftval / 2

    T result;

    if (xi_lr * (left + right) < left )
    {
        //sample from left triangle
        result = 1 - sample_linear(xi_linear);
    }
    else
    {
        // sample from right triangle
        result = sample_linear(xi_linear);
    }

    ENSURE(0 <= result && result <= 1);
    return result;
}

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
template<class FloatType>
__device__  inline FloatType * sample_dcdf(
        FloatType   xi,
        FloatType  *iter,
        FloatType  *last)
{
    REQUIRE(iter-last > 0);
    REQUIRE(xi >= 0 && xi < 1);

    // Do a binary search on the CDF
    iter = cuda::utility::lower_bound(iter, last, xi);

    // Return the value
    ENSURE(iter != last);
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
template<class FloatType>
__device__  inline FloatType * sample_small_dcdf(
        FloatType   xi,
        FloatType  *iter,
        FloatType  *last)
{
    REQUIRE(iter - last > 0);
    REQUIRE(xi >= 0 && xi < 1);

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
template<class FloatType>
__device__  inline FloatType * sample_smallrev_dcdf(
        FloatType   xi,
        FloatType  *first,
        FloatType  *iter)
{
    REQUIRE(first - iter > 0);
    REQUIRE(xi >= 0 && xi < 1);

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
__device__  inline T sample_linear(T xi)
{
    REQUIRE(0 <= xi && xi < 1);

    return sqrt(xi);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a linear distribution on [0,1) with arbitrary y values
 *
 * This pulls two random numbers from the generator and passes them to the
 * underlying sampler function.
 */
template<class RNG_State>
__device__  inline float sample_linear(RNG_State  *rng,
                                       const float left,
                                       const float right)
{
    return sample_linear(
            curand_uniform(rng), curand_uniform(rng),
            left, right);
}

template<class RNG_State>
__device__  inline double sample_linear(RNG_State *rng, const double left,
                                        const double right)
{
    return sample_linear(
            curand_uniform_double(rng), curand_uniform_double(rng),
            left, right);
}

//---------------------------------------------------------------------------//
//!\brief Sample an isotropic angular distribution
__device__ inline void sample_isotropic(cuda_utils::Space_Vector &omega,
                                        RNG_State_t              *rng)
{
    omega.z         = 1.0 - 2.0 * curand_uniform_double(rng);
    double phi      = cuda::constants::two_pi *
                      curand_uniform_double(rng);
    double sintheta = sqrt(1.0 - omega.z * omega.z);

    omega.x = sintheta * cos(phi);
    omega.y = sintheta * sin(phi);
}

//---------------------------------------------------------------------------//

} // end namespace cuda_mc::sampler

} // end namespace cuda_mc

#endif // cuda_mc_Sampler_i_cuh

//---------------------------------------------------------------------------//
//                 end of Sampler.i.cuh
//---------------------------------------------------------------------------//
