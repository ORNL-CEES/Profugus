//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Vector_Functions.hh
 * \author Thomas M. Evans
 * \date   Tuesday April 29 14:18:9 2014
 * \brief  General vector functions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Vector_Functions_hh
#define utils_Vector_Functions_hh

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "Definitions.hh"
#include "Vector_Lite.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// SPACE_VECTOR FUNCTIONS
//---------------------------------------------------------------------------//

// Transform a space-vector through (theta, phi) in cartesian space.
void cartesian_vector_transform(double costheta, double phi,
                                def::Space_Vector &vector);

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the magnitude of a space-vector.
 *
 * The magnitude of a vector is defined,
 * \f[
   |\mathbf{V}| = \sqrt{V_x^2 + V_y^2 + V_z^2}
 * \f]
 *
 * \return vector magnitude
 */
inline double vector_magnitude(const def::Space_Vector &vector)
{
    return std::sqrt(vector[def::X] * vector[def::X] +
                     vector[def::Y] * vector[def::Y] +
                     vector[def::Z] * vector[def::Z]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the dot product of two space-vectors.
 *
 * The dot product is defined as
 * \f[
 *  \mathbf{v} \cdot \mathbf{w} = v_x*w_x + v_y*w_y + v_z*w_z \, .
 * \f]
 *
 * \return dot product
 */
inline double dot_product(const def::Space_Vector &v,
                          const def::Space_Vector &w)
{
    return v[def::X]*w[def::X] + v[def::Y]*w[def::Y] + v[def::Z]*w[def::Z];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Normalize a vector so that its magnitude is one.
 *
 * Convert a vector into a unit-vector using,
 * \f[
   \hat{\mathbf{V}} = \frac{\mathbf{V}}{|\mathbf{V}|},
 * \f]
 *
 * This function is unrolled for efficiency, which is why we don't have a
 * general normalize function.
 */
inline void vector_normalize(def::Space_Vector &vector)
{
    register double norm = 1.0 / vector_magnitude(vector);
    vector[def::X] *= norm;
    vector[def::Y] *= norm;
    vector[def::Z] *= norm;

    Ensure (soft_equiv(vector_magnitude(vector), 1.0, 1.0e-6));
}

} // end namespace profugus

#endif // utils_Vector_Functions_hh

//---------------------------------------------------------------------------//
//              end of utils/Vector_Functions.hh
//---------------------------------------------------------------------------//
