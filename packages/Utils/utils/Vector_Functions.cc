//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Vector_Functions.cc
 * \author Thomas M. Evans
 * \date   Mon Aug 22 11:21:49 2011
 * \brief  Vector function definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include "Vector_Functions.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Transform a space-vector from \f$\Omega\rightarrow\Omega'\f$ through
 * \f$(\theta, \phi)\f$ in cartesian space.
 *
 * The direction cosines are given by the formula:
 * \f[
 * \mathbf{\Omega} = \Omega_{x}\mathbf{i} + \Omega_{y}\mathbf{j}
 * + \Omega_{z}\mathbf{k}
 * \f]
 * When scattered through an angle \f$(\theta,\phi)\f$, the new direction
 * cosines are:
 * \f[
 * \Omega_{x}' = \Omega_{x}\cos\theta + \Omega_{x}\Omega_{z}\sin\theta
 * \cos\phi / \alpha - \Omega_{y}\sin\theta\sin\phi/\alpha
 * \f]
 * \f[
 * \Omega_{y}' = \Omega_{y}\cos\theta + \Omega_{y}\Omega_{z}\sin\theta
 * \cos\phi / \alpha + \Omega_{x}\sin\theta\sin\phi/\alpha
 * \f]
 * \f[
 * \Omega_{z}' = \Omega_{z}\cos\theta - \alpha\sin\theta\cos\phi
 * \f]
 * where
 * \f[
 * \alpha = \sqrt{1-\Omega_{z}^{2}}
 * \f]
 *
 * \param costheta cosine of polar scattering angle
 * \param phi azimuthal scattering angle
 * \param vector on input, the initial direction; on output the final
 * direction transformed through \f$(\theta,\phi)\f$
 *
 * \pre the vector must be a unit-vector (magnitude == 1)
 */
void cartesian_vector_transform(double             costheta,
                                double             phi,
                                def::Space_Vector &vector)
{
    using def::X; using def::Y; using def::Z;

    Require (soft_equiv(vector_magnitude(vector), 1.0, 1.0e-6));

    // cos/sin factors
    const double cosphi   = std::cos(phi);
    const double sinphi   = std::sin(phi);
    const double sintheta = std::sqrt(1.0 - costheta * costheta);

    // make a copy of the old direction
    def::Space_Vector old = vector;

    // calculate alpha
    const double alpha = std::sqrt(1.0 - old[Z] * old[Z]);

    // now transform into new cooordinate direction; degenerate case first
    if (alpha < 1.e-6)
    {
        vector[X] = sintheta * cosphi;
        vector[Y] = sintheta * sinphi;
        vector[Z] = (old[Z] < 0.0 ? -1.0 : 1.0) * costheta;
    }

    // do standard transformation
    else
    {
        // calculate inverse of alpha
        register const double inv_alpha = 1.0 / alpha;

        // calculate new z-direction
        vector[Z] = old[Z] * costheta - alpha * sintheta * cosphi;

        // calculate new x-direction
        vector[X] = old[X] * costheta + inv_alpha * (
            old[X] * old[Z] * sintheta * cosphi - old[Y] * sintheta * sinphi);

        // calculate new y-direction
        vector[Y] = old[Y] * costheta + inv_alpha * (
            old[Y] * old[Z] * sintheta * cosphi + old[X] * sintheta * sinphi);
    }

    // normalize the particle to avoid roundoff errors
    vector_normalize(vector);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Vector_Functions.cc
//---------------------------------------------------------------------------//
