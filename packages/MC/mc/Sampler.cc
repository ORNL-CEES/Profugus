//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Sampler.cc
 * \author Thomas M. Evans
 * \date   Friday May 2 10:28:14 2014
 * \brief  Sampling function definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Sampler.hh"
#include "rng/RNG.hh"

namespace profugus
{

namespace sampler
{

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a discrete CDF.
 *
 * Given a CDF defined over the range \f$[0,N_b)\f$, where \f$N_b\f$ is the
 * number of bins, sample and return the bin of the sampled value. The domain
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
int sample_discrete_CDF(int      nb,
                        const T *c,
                        const T  ran)
{
    REQUIRE(nb > 0);
    REQUIRE(soft_equiv(static_cast<double>(c[nb - 1]), 1.0, 1.0e-6));
    REQUIRE(ran >= 0 && ran <= 1);

    // do a binary search on the CDF
    const T *ptr = std::lower_bound(c, c + nb, ran);

    // return the value
    ENSURE(ptr - c >= 0 && ptr - c < nb);
    return ptr - c;
}

// Explicit instantiations on float and double
template int sample_discrete_CDF(int nb, const float *c , const float  ran);
template int sample_discrete_CDF(int nb, const double *c, const double ran);

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
template<class RNG>
double sample_watt(RNG   &rng,
                   double a,
                   double b)
{
    using std::sqrt;
    using std::log;

    REQUIRE(a > 0);
    REQUIRE(b > 0);

    // the algorithm implemented here is the same as in MCNP5

    const double t1 = 1 + 0.125 * a * b;
    const double c1 = t1 + sqrt(t1 * t1 - 1) - 1;
    const double c2 = a * (c1 + 1);
    const double c3 = b * c2;

    while (1)
    {
        double r1 = -log(rng.ran());
        double r2 = -log(rng.ran());
        double c4 = r2 - c1 * (1 + r1);
        if (c4 * c4 <= c3 * r1)
            return 1.0e6 * c2 * r1;  // in eV
    }
}

//---------------------------------------------------------------------------//

// Explicit instantiation for SPRNG
template double sample_watt(RNG &rng, double a, double b);

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
T sample_linear(
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
        result = 1 - sample_linear<T>(xi_linear);
    }
    else
    {
        // sample from right triangle
        result = sample_linear<T>(xi_linear);
    }

    ENSURE(0 <= result && result <= 1);
    return result;
}

// Explicit instantiation
template double sample_linear<double>(double, double, double, double);
template float  sample_linear<float >(float, float, float, float);

//---------------------------------------------------------------------------//
} // end namespace sampler

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Sampler.cc
//---------------------------------------------------------------------------//
