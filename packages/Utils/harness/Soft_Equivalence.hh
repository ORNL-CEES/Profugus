//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/harness/Soft_Equivalence.hh
 * \author Thomas M. Evans and Seth R Johnson
 * \date   Wed Jan  2 11:56:34 2008
 * \brief  Soft_Equivalence functions for floating point comparisons.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_harness_Soft_Equivalence_hh
#define Utils_harness_Soft_Equivalence_hh

#include <cmath>
#include "DBC.hh"

namespace profugus
{
namespace detail
{

//---------------------------------------------------------------------------//
/*!
 * \struct softeq_traits
 * \brief Provide relative errors for soft_equiv based on type.
 *
 * This also gives compile-time checking for bad values.
 */
template<typename T>
struct softeq_traits
{
    typedef T value_type;
    //! Default relative error
    static constexpr value_type rel_prec()
    {
        static_assert(sizeof(T) == 0, "Invalid type for softeq!");
        return T();
    }

    //! Default absolute error
    static constexpr value_type abs_thresh()
    {
        static_assert(sizeof(T) == 0, "Invalid type for softeq!");
        return T();
    }
};

template<>
struct softeq_traits<double>
{
    typedef double value_type;
    static constexpr value_type rel_prec() { return 1.0e-12; }
    static constexpr value_type abs_thresh() { return 1.0e-14; }
};

template<>
struct softeq_traits<float>
{
    typedef float value_type;
    static constexpr value_type rel_prec() { return 1.0e-6f; }
    static constexpr value_type abs_thresh() { return 1.0e-8f; }
};

//---------------------------------------------------------------------------//
/*!
 * \struct precision_type
 * \brief Get a "least common denominator" for soft comparisons.
 */
template<typename Expected_T, typename Actual_T>
struct precision_type
{
    // Not available
    typedef int type;
};

template<>
struct precision_type<double,float>
{
    typedef float type;
};

template<>
struct precision_type<float,float>
{
    typedef float type;
};

template<>
struct precision_type<float,double>
{
    typedef float type;
};

template<>
struct precision_type<double,double>
{
    typedef double type;
};

// In the case where user gives literal 0 as the reference value
template<>
struct precision_type<int,double>
{
    typedef double type;
};

}

//===========================================================================//
// SCALAR SOFT EQUIVALENCE FUNCTIONS
//===========================================================================//
/*!
 * \brief Compare two floating point scalars for equivalence to a specified
 * tolerance.
 *
 * \param value scalar floating point value
 *
 * \param reference scalar floating point reference to which value is
 * compared
 *
 * \param rel_precision tolerance of relative error (default 1.0e-12)
 *
 * \param abs_threshold threshold for absolute error when comparing to zero
 *                      (default 1.0e-14)
 *
 * \return true if values are the same within relative error specified by
 * precision, false if otherwise
 */
template<typename Expected_T, typename Actual_T>
inline bool soft_equiv(
    const Actual_T   value,
    const Expected_T reference,
    const typename detail::precision_type<
                                Expected_T,Actual_T>::type rel_precision,
    const typename detail::precision_type<
                                Expected_T,Actual_T>::type abs_threshold)
{
    using std::fabs;

    typedef typename detail::precision_type<Expected_T,Actual_T>::type
        result_type;

    // Promote/demote value type
    result_type val = static_cast<result_type>(value);
    result_type ref = static_cast<result_type>(reference);

    // Typical case: relative error comparison to reference
    if (fabs(val - ref) < rel_precision * fabs(ref))
    {
        return true;
    }

    // If one is within the absolute threshold of zero, and the other within
    // relative of zero, they're equal
    if ((fabs(ref) < abs_threshold) && (fabs(val) < rel_precision))
    {
        return true;
    }
    if ((fabs(val) < abs_threshold) && (fabs(ref) < rel_precision))
    {
        return true;
    }

    return false;
}

//! Soft equiv with default absolute precision
template<typename Expected_T, typename Actual_T>
inline bool soft_equiv(
        const Actual_T   value,
        const Expected_T reference,
        const typename detail::precision_type<Expected_T,Actual_T>::type rel)
{
    typedef typename detail::precision_type<Expected_T,Actual_T>::type prec_t;
    return soft_equiv(value, reference,
                      rel,
                      detail::softeq_traits<prec_t>::abs_thresh());
}

//! Soft equiv with default absolute and relative precision
template<typename Expected_T, typename Actual_T>
inline bool soft_equiv(
        const Actual_T   value,
        const Expected_T reference)
{
    typedef typename detail::precision_type<Expected_T,Actual_T>::type prec_t;
    return soft_equiv(value, reference,
                      detail::softeq_traits<prec_t>::rel_prec(),
                      detail::softeq_traits<prec_t>::abs_thresh());
}

} // end namespace profugus

#endif // Utils_harness_Soft_Equivalence_hh

//---------------------------------------------------------------------------//
//              end of harness/Soft_Equivalence.hh
//---------------------------------------------------------------------------//
