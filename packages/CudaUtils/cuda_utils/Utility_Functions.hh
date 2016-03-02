//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Utility_Functions.hh
 * \author Steven Hamilton
 * \date   Tue Dec 15 16:25:39 2015
 * \brief  Utility_Functions class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Utility_Functions_hh
#define CudaUtils_cuda_utils_Utility_Functions_hh

#include "Definitions.hh"

namespace cuda
{

namespace utility
{

// On-device lower bound binary search
template <class T>
__host__ __device__ inline const T * lower_bound(const T *first,
                                                 const T *last,
                                                 const T &val )
{
    const T *it;
    int count, step;
    count = last - first;
    while( count > 0 )
    {
        step = count / 2;
        it = first + step;
        if( *it < val )
        {
            first = ++it;
            count -= step+1;
        }
        else
        {
            count = step;
        }
    }
    return first;
}

__host__ __device__ inline bool soft_equiv( double a,
                                            double b,
                                            double tol = 1.0e-12 )
{
    if( abs(a - b) < tol )
        return true;

    if( abs(a) < tol && abs(b) < tol )
        return true;

    return false;
}


// Cuda doesn't natively support atomic operations on doubles
double __device__ inline atomic_add_double( volatile double * const dest,
                                            const double val )
{
    union U
    {
        unsigned long long int i;
        double t;
    } assume, oldval, newval;

    oldval.t = *dest;

    do
    {
        assume.i = oldval.i;
        newval.t = assume.t + val;
        oldval.i = atomicCAS( (unsigned long long int*)dest,
                              assume.i, newval.i );
    } while (assume.i != oldval.i);

    return oldval.t;
}


//---------------------------------------------------------------------------//
// SPACE_VECTOR FUNCTIONS
//---------------------------------------------------------------------------//

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
__device__ inline double vector_magnitude(const cuda::Space_Vector &vector)
{
    return sqrt(vector.x * vector.x +
                vector.y * vector.y +
                vector.z * vector.z);
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
__device__ inline double dot_product(const cuda::Space_Vector &v,
                                     const cuda::Space_Vector &w)
{
    return v.x*w.x + v.y*w.y + v.z*w.z;
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
__device__ inline void vector_normalize(cuda::Space_Vector &vector)
{
    double norm     = 1.0 / vector_magnitude(vector);
    vector.x *= norm;
    vector.y *= norm;
    vector.z *= norm;

    ENSURE(soft_equiv(vector_magnitude(vector), 1.0, 1.0e-6));
}

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
__device__ inline void cartesian_vector_transform(double costheta, double phi,
                                                  cuda::Space_Vector &vector)
{
    REQUIRE(soft_equiv(vector_magnitude(vector), 1.0, 1.0e-6));

    // cos/sin factors
    const double cosphi   = cos(phi);
    const double sinphi   = sin(phi);
    const double sintheta = sqrt(1.0 - costheta * costheta);

    // make a copy of the old direction
    cuda::Space_Vector old = vector;

    // calculate alpha
    const double alpha = sqrt(1.0 - old.z * old.z);

    // now transform into new cooordinate direction; degenerate case first
    if (alpha < 1.e-6)
    {
        vector.x = sintheta * cosphi;
        vector.y = sintheta * sinphi;
        vector.z = (old.z < 0.0 ? -1.0 : 1.0) * costheta;
    }

    // do standard transformation
    else
    {
        // calculate inverse of alpha
        const double inv_alpha = 1.0 / alpha;

        // calculate new z-direction
        vector.z = old.z * costheta - alpha * sintheta * cosphi;

        // calculate new x-direction
        vector.x = old.x * costheta + inv_alpha * (
            old.x * old.z * sintheta * cosphi - old.y * sintheta * sinphi);

        // calculate new y-direction
        vector.y = old.y * costheta + inv_alpha * (
            old.y * old.z * sintheta * cosphi + old.x * sintheta * sinphi);
    }

    // normalize the particle to avoid roundoff errors
    vector_normalize(vector);

}




} // end namespace utility

} // end namespace cuda

#endif // CudaUtils_cuda_utils_Utility_Functions_hh

//---------------------------------------------------------------------------//
// end of CudaUtils/cuda_utils/Utility_Functions.hh
//---------------------------------------------------------------------------//
