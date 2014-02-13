//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Container_Functions.hh
 * \author Gregory G. Davidson
 * \date   Thu Jan 09 11:56:28 2014
 * \brief  Various container-related functions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Container_Functions_hh
#define utils_Container_Functions_hh

#include "Soft_Comparators.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// REMOVE SUBSETS
//---------------------------------------------------------------------------//
// Return container 1 with the elements in container 2 removed (if they
// exist).  Containers do not need to be sorted.
template<class Container1, class Container2>
void remove_subset(Container1       &cont1,
                   const Container2 &cont2);

//---------------------------------------------------------------------------//
// Return container 1 with the elements in container 2 removed, using Comp as
// the comparator.  Containers do not need to be sorted.
template<class Container1, class Container2, class BinaryOp>
void remove_subset(Container1       &cont1,
                   const Container2 &cont2,
                   BinaryOp          comp);

//---------------------------------------------------------------------------//
// CONTAINER MERGE FUNCTIONS
//---------------------------------------------------------------------------//
// Merge the elements in container 2 into container 1.  Containers must be
// sorted
template<class Container1, class Container2>
void container_merge(Container1       &cont_1,
                     const Container2 &cont_2);

// Merge the elements in container 2 into container 1, within epsilon.
// Containers must be sorted.
template<class Container1, class Container2>
void container_merge(Container1                            &cont1,
                     const Container2                      &cont2,
                     const typename Container1::value_type &eps);

// Merges the elements in container 2 into container 1, using operator
// EqlComp.  Containers must have been sorted with LessComp.
template<class Container1, class Container2, class EqualOp, class LessOp>
void container_merge(Container1       &cont1,
                     const Container2 &cont2,
                     EqualOp           eql_comp,
                     LessOp            less_comp);

//---------------------------------------------------------------------------//
// ASSOCIATIVE CONTAINER FUNCTIONS
//---------------------------------------------------------------------------//

namespace CopyType
{

//! An enumeration type indicating whether to copy the keys or the values.
enum CT{ KEY, VALUE };

} // end namespace profugus::CopyType

// Copying from an associative container to a single-value container
template<class AC_ConstIterator, class SV_Iterator>
void copy(AC_ConstIterator ac_begin,
          AC_ConstIterator ac_end,
          SV_Iterator      sv_begin,
          CopyType::CT     copy_type = CopyType::VALUE);

//---------------------------------------------------------------------------//
// Fill a map from two provided single-value containers
template<class SV_Iterator1, class SV_Iterator2, class AS_Container>
void fill_map(SV_Iterator1  key_iter_begin,
              SV_Iterator1  key_iter_end,
              SV_Iterator2  value_iter_begin,
              AS_Container &assoc_cont);

//---------------------------------------------------------------------------//
// Create a multimap from two provided single-value containers
template<class SV_Iterator1, class SV_Iterator2, class AS_Container>
void fill_multimap(SV_Iterator1  key_iter_begin,
                   SV_Iterator1  key_iter_end,
                   SV_Iterator2  value_iter_begin,
                   AS_Container &assoc_cont);

//---------------------------------------------------------------------------//
// Split an associative container into two single-value containers holding the
// keys and values in order
template<class AC_Iterator, class SV_Iterator1, class SV_Iterator2>
void split_map(AC_Iterator  map_iter_begin,
               AC_Iterator  map_iter_end,
               SV_Iterator1 key_iter,
               SV_Iterator2 value_iter);

//---------------------------------------------------------------------------//
// TRIM FUNCTIONS
//---------------------------------------------------------------------------//
// Trim a sorted container by cutting off everything that lies outside the
// given bounds
template<class Container>
void trim(Container                            &container,
          const typename Container::value_type &lower_bound,
          const typename Container::value_type &upper_bound);

//---------------------------------------------------------------------------//
// Trim a sorted container by cutting off everything that lies outside the
// given bounds within epsilon
template<class Container>
void trim(Container                            &container,
          const typename Container::value_type &lower_bound,
          const typename Container::value_type &upper_bound,
          const typename Container::value_type &eps);

//---------------------------------------------------------------------------//
// Trim a container using EqlComp that has been sorted with LessComp.
template<class Container, class EqualOp, class LessOp>
void trim(Container                            &container,
          const typename Container::value_type &lower_bound,
          const typename Container::value_type &upper_bound,
          EqualOp                               equal_comp,
          LessOp                                less_comp);

//---------------------------------------------------------------------------//
// UNIQUE-ING FUNCTIONS
//---------------------------------------------------------------------------//
// Make the elements unique.  Does not assume the container is sorted.
template<class Container>
void make_unique(Container &container);

//---------------------------------------------------------------------------//
// Make the elements unique within epsilon.
template<class Container>
void make_unique(Container                            &container,
                 const typename Container::value_type &epsilon);

//---------------------------------------------------------------------------//
// Make the elements unique given equality and sorting operators.
template<class Container, class EqualOp, class LessOp>
void make_unique(Container &container,
                 EqualOp    equal_comp,
                 LessOp     less_comp);

//---------------------------------------------------------------------------//
// NUMERICAL FUNCTIONS
//---------------------------------------------------------------------------//
// Return the arithmetic mean of the elements
template<class InputIterator>
typename std::iterator_traits<InputIterator>::value_type
mean(InputIterator begin,
     InputIterator end);

//---------------------------------------------------------------------------//
// Return the median of the elements in a container
template<class InputIterator>
typename std::iterator_traits<InputIterator>::value_type
median(InputIterator begin,
       InputIterator end);

//---------------------------------------------------------------------------//
// Truncate a vector to [beg, end), renormalizing it
template<class Container>
typename Container::value_type
truncate_and_norm(typename Container::iterator  begin,
                  typename Container::iterator  end,
                  Container                    &t);

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE FUNCTION DEFINITIONS
//---------------------------------------------------------------------------//

#include "Container_Functions.i.hh"

#endif // utils_Container_Functions_hh

//---------------------------------------------------------------------------//
//              end of utils/Container_Functions.hh
//---------------------------------------------------------------------------//
