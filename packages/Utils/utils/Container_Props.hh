//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Container_Props.hh
 * \author Gregory G. Davidson
 * \date   Thu Jan 09 11:13:13 2014
 * \brief  Container_Props class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Container_Props_hh
#define utils_Container_Props_hh

#include "Soft_Comparators.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// IS_SORTED FUNCTIONS
//---------------------------------------------------------------------------//
// Return whether a container is sorted.
template<class InputIterator>
bool is_sorted(InputIterator begin,
               InputIterator end);

//---------------------------------------------------------------------------//
// Return whether a container is sorted according to the given comparator.
template<class InputIterator, class BinaryOp>
bool is_sorted(InputIterator begin,
               InputIterator end,
               BinaryOp      comp);

//---------------------------------------------------------------------------//
// IS_SUBSET FUNCTIONS
//---------------------------------------------------------------------------//
// Return whether cont2 is a subset of cont1
template<class Container1, class Container2>
bool is_subset(const Container1 &cont1,
               const Container2 &cont2);

// Return whether cont2 is a subset of cont1 within the given epsilon
template<class Container1, class Container2>
bool is_subset(const Container1                      &cont1,
               const Container2                      &cont2,
               const typename Container1::value_type &eps);

//---------------------------------------------------------------------------//
// IS_UNIQUE FUNCTIONS
//---------------------------------------------------------------------------//
// Return whether all elements in a container are unique
template<class InputIterator>
bool is_unique(InputIterator begin,
               InputIterator end);

//---------------------------------------------------------------------------//
// Return whether all elements in a container are unique within epsilon.
template<class InputIterator>
bool is_unique(InputIterator                          begin,
               InputIterator                          end,
               const typename ITT<InputIterator>::vt &epsilon);

//---------------------------------------------------------------------------//
// Return whether all elements in a container are unique given equality and
// sorting operators
template<class InputIterator, class EqualOp, class LessOp>
bool is_unique(InputIterator begin,
               InputIterator end,
               EqualOp       equal_comp,
               LessOp        less_comp);

//---------------------------------------------------------------------------//
// NUMERICAL FUNCTIONS
//---------------------------------------------------------------------------//
// Check a range of elements to ensure its contents are all positive (or zero)
template<class InputIterator>
bool is_positive(InputIterator begin,
                 InputIterator end);

//---------------------------------------------------------------------------//
// Check a container to ensure its contents are negative (or zero).
template<class InputIterator>
bool is_negative(InputIterator begin,
                 InputIterator end);

//---------------------------------------------------------------------------//
// Check a container to ensure its contents are non-negative.
template<class InputIterator>
bool is_non_negative(InputIterator begin,
                     InputIterator end);

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE FUNCTION DEFINITIONS
//---------------------------------------------------------------------------//

#include "Container_Props.i.hh"

#endif // utils_Container_Props_hh

//---------------------------------------------------------------------------//
//                      end of utils/Container_Props.hh
//---------------------------------------------------------------------------//
