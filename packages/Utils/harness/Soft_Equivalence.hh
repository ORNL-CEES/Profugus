//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/Soft_Equivalence.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 11:56:34 2008
 * \brief  Soft_Equivalence functions for floating point comparisons.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef harness_Soft_Equivalence_hh
#define harness_Soft_Equivalence_hh

#include <cmath>
#include <iterator>
#include "DBC.hh"

namespace profugus
{

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
 * \param precision tolerance of relative error (default 1.0e-12)
 *
 * \return true if values are the same within relative error specified by
 * precision, false if otherwise
 *
 * \todo Should we be using numeric_limits instead of hard coded vales for
 *       e-12 and e-14?
 */
template<class FPT>
inline bool soft_equiv(const FPT &value,
                       const FPT &reference,
                       const FPT  precision = 1.0e-12)
{
    using std::fabs;
    bool passed = false;

    if (fabs(value - reference) < precision * fabs(reference))
            passed = true;
    else
            passed = false;

    // second chance for passing if reference is within machine error of zero
    if (!passed && (fabs(reference) < 1.0e-14))
            if (fabs(value) < precision)
            passed = true;

    // third chance for passing if value is within machine error of zero
    if (!passed && (fabs(value) < 1.0e-14))
        if (fabs(reference) < precision)
            passed = true;

    return passed;
}

//---------------------------------------------------------------------------//

template<>
inline bool soft_equiv(const int &value,
                       const int &reference,
                       const int  precision)
{
    Insist (0, "Can't do a soft compare with integers!");
    return false;
}

//===========================================================================//
// FIELD SOFT EQUIVALENCE FUNCTIONS
//===========================================================================//
/*!
 * \brief Compare two floating point fields for equivalence to a specified
 * tolerance.
 *
 * \param value  floating point field of values
 * \param value_end End of floating point field
 * \param ref floating point field to which values are compared
 * \param ref_end End of floating point reference field
 *
 * \param precision tolerance of relative error (default 1.0e-12)
 *
 * \return true if values are the same within relative error specified by
 * precision and the fields are the same size, false if otherwise
 *
 * The field soft_equiv check is an element-by-element check of two
 * single-dimension fields.  The precision is the same type as the value
 * field.  The value and reference fields must have STL-type iterators.  The
 * value-types of both fields must be the same or a compile-time error will
 * result.
 */
template<class Value_Iterator, class Ref_Iterator>
inline bool soft_equiv(
    Value_Iterator value,
    Value_Iterator value_end,
    Ref_Iterator   ref,
    Ref_Iterator   ref_end,
    const typename std::iterator_traits<Value_Iterator>::value_type precision
    = 1.0e-12)
{
    using std::distance;

    bool passed = true;

    // first check that the sizes are equivalent
    if (distance(value, value_end) != distance(ref, ref_end))
    {
        passed = false;
    }

    // if the sizes are the same, loop through and check each element
    else
    {
        while (value != value_end && passed == true)
        {
            passed = soft_equiv(*value, *ref, precision);
            value++;
            ref++;
        }
    }

    return passed;
}

} // end namespace profugus

#endif // harness_Soft_Equivalence_hh

//---------------------------------------------------------------------------//
//              end of harness/Soft_Equivalence.hh
//---------------------------------------------------------------------------//
