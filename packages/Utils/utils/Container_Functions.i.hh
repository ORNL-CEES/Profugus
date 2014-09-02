//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Container_Functions.i.hh
 * \author Gregory G. Davidson
 * \date   Thu Jan 09 11:56:28 2014
 * \brief  Inline function definitions of of container functions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Container_Functions_i_hh
#define utils_Container_Functions_i_hh

#include <numeric>

#include "harness/Soft_Equivalence.hh"
#include "Container_Props.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// SUBSET FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Remove from \a cont1 the elements in \a cont2.
 */
template<class Container1, class Container2>
void remove_subset(Container1       &cont1,
                   const Container2 &cont2)
{
    // Create the comparator
    typedef typename Container1::value_type             val_type;
    typedef typename Compare_Traits<val_type>::eql_comp eql_comp;

    // Remove the subset
    remove_subset(cont1, cont2, eql_comp());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Remove from \a cont1 the elements in \a cont2, using \a comp as the
 *        comparator.
 */
template<class Container1, class Container2, class BinaryOp>
void remove_subset(Container1       &cont1,
                   const Container2 &cont2,
                   BinaryOp          comp)
{
    // Loop over container 2
    for(typename Container2::const_iterator iter = cont2.begin(),
                                        iter_end = cont2.end();
        iter != iter_end; ++iter)
    {
        cont1.erase( std::remove_if(cont1.begin(), cont1.end(),
                                    bind2nd(comp, *iter)), cont1.end() );
    }
}

//---------------------------------------------------------------------------//
// CONTAINER MERGE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Merge \a cont2 into \a cont1.
 *
 * \pre Containers \a cont_1 and \a cont_2 must be sorted.
 */
template<class Container1, class Container2>
void container_merge(Container1       &cont_1,
                     const Container2 &cont_2)
{
    typedef typename Container1::value_type              val_type;
    typedef typename Compare_Traits<val_type>::less_comp less_comp;
    typedef typename Compare_Traits<val_type>::eql_comp  eql_comp;

    REQUIRE(profugus::is_sorted(cont_1.begin(), cont_1.end()));
    REQUIRE(profugus::is_sorted(cont_2.begin(), cont_2.end()));

    container_merge(cont_1, cont_2, eql_comp(), less_comp());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Merge the elements in \a cont2 into \a cont1 within \a epsilon.
 *
 * \pre Containers \a cont_1 and \a cont_2 must be sorted within epsilon.
 */
template<class Container1, class Container2>
void container_merge(Container1                            &cont_1,
                     const Container2                      &cont_2,
                     const typename Container1::value_type &eps)
{
    typedef typename Container1::value_type  value_type;

    REQUIRE( profugus::is_sorted(cont_1.begin(), cont_1.end(),
                                  SoftIsLessEqual<value_type>(eps)) );
    REQUIRE( profugus::is_sorted(cont_2.begin(), cont_2.end(),
                                  SoftIsLessEqual<value_type>(eps)) );

    container_merge(cont_1, cont_2, SoftIsEqual<value_type>(eps),
                    SoftIsLessEqual<value_type>(eps));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Merge elements in \a cont2 into \a cont1, using operator \a
 *        eql_comp.
 *
 * Example 1: Union of {0.0, 2.0, 3.0, 4.0} and {0.0, 4.0} =
 *                                                         {0.0, 2.0, 3.0, 4.0}
 * Example 2: Union of {0.0, 1.0, 4.0} and {2.0, 4.0}      =
 *                                                         {0.0, 1.0, 2.0, 4.0}
 *
 * \note Containers must be sorted with operator \a less_comp.
 */
template<class Container1, class Container2, class EqualOp, class LessOp>
void container_merge(Container1       &cont_1,
                     const Container2 &cont_2,
                     EqualOp           eql_comp,
                     LessOp            less_comp)
{
    REQUIRE( profugus::is_sorted(cont_1.begin(), cont_1.end(), less_comp) );
    REQUIRE( profugus::is_sorted(cont_2.begin(), cont_2.end(), less_comp) );

    // Iterators into the two containers to be merged
    typename Container1::iterator iter_1       = cont_1.begin();
    typename Container2::const_iterator iter_2 = cont_2.begin();

    // Loop until both iterators have advanced to the end
    while(iter_1 != cont_1.end() || iter_2 != cont_2.end() )
    {
        // If iter_1 is at cont_1.end(), insert remainder of cont_2 and
        // advance iter_2 and iter_1 to end
        if(iter_1 == cont_1.end())
        {
            cont_1.insert(cont_1.end(), iter_2, cont_2.end());
            iter_2 = cont_2.end();
            iter_1 = cont_1.end();
        }
        // If iter_2 is at cont_2.end(), advance iter_1 to end
        else if(iter_2 == cont_2.end())
        {
            iter_1 = cont_1.end();
        }
        // If iter_1 and iter_2 point to a value that is "equal",
        // advance both iterators
        else if( eql_comp(*iter_1, *iter_2) )
        {
            ++iter_1, ++iter_2;
        }
        // If iter_1 points to a value less than iter_2, advance iter_1
        else if( less_comp(*iter_1, *iter_2) )
        {
            ++iter_1;
        }
        // If we get down here, *iter_2 must be smaller than *iter_1.  Insert
        // *iter_2 at position iter_1.  Reset iter_1 and advance iter_2
        else
        {
            CHECK( less_comp(*iter_2, *iter_1) );
            iter_1 = cont_1.insert(iter_1, *iter_2);
            ++iter_2;
        }
    }
}

//---------------------------------------------------------------------------//
// ASSOCIATIVE CONTAINER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Copies the key or value or an associative container into a
 *        single-value container.
 *
 * The resultant single-value container will be sorted in the same order as
 * the associative container. It must be the same size (or larger) as the
 * source range.
 */
template<class AC_ConstIterator, class SV_Iterator>
void copy(AC_ConstIterator ac_begin,
          AC_ConstIterator ac_end,
          SV_Iterator      sv_begin,
          CopyType::CT     copy_type)
{
    if(copy_type == CopyType::KEY)
    {
        while(ac_begin != ac_end)
        {
            *sv_begin = ac_begin->first;
            ++sv_begin, ++ac_begin;
        }
    }
    else
    {
        while(ac_begin != ac_end)
        {
            *sv_begin = ac_begin->second;
            ++sv_begin, ++ac_begin;
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Fill a map from two provided single-value containers.
 */
template<class SV_Iterator1, class SV_Iterator2, class AS_Container>
void fill_map(SV_Iterator1  key_iter_begin,
              SV_Iterator1  key_iter_end,
              SV_Iterator2  value_iter_begin,
              AS_Container &assoc_cont)
{
    REQUIRE(is_unique(key_iter_begin, key_iter_end));
    // Make sure none of the keys are already in the map
#ifdef REQUIRE_ON
    for(SV_Iterator1 iter = key_iter_begin; iter != key_iter_end; ++iter)
    {
        REQUIRE(assoc_cont.count(*iter) == 0);
    }
#endif

    // Loop over the keys
    while(key_iter_begin != key_iter_end)
    {
        CHECK(assoc_cont.count(*key_iter_begin) == 0);

        // Insert the value into the associative container
        assoc_cont.insert( std::make_pair(*key_iter_begin, *value_iter_begin) );

        // Update the iterator
        ++key_iter_begin, ++value_iter_begin;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a multimap from two provided single-value containers.
 */
template<class SV_Iterator1, class SV_Iterator2, class AS_Container>
void fill_multimap(SV_Iterator1  key_iter_begin,
                   SV_Iterator1  key_iter_end,
                   SV_Iterator2  value_iter_begin,
                   AS_Container &assoc_cont)
{
    // Loop over the keys
    while(key_iter_begin != key_iter_end)
    {
        // Insert the value into the associative container
        assoc_cont.insert( std::make_pair(*key_iter_begin, *value_iter_begin) );

        // Update the iterator
        ++key_iter_begin, ++value_iter_begin;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Split a map into two single-value containers holding the keys and
 *        values of the map.
 *
 * The keys and values containers must be large enough to hold the source
 * range.
 */
template<class AC_Iterator, class SV_Iterator1, class SV_Iterator2>
void split_map(AC_Iterator  map_iter_begin,
               AC_Iterator  map_iter_end,
               SV_Iterator1 key_iter,
               SV_Iterator2 value_iter)
{
    while(map_iter_begin != map_iter_end)
    {
        *key_iter   = map_iter_begin->first;
        *value_iter = map_iter_begin->second;

        ++map_iter_begin;
        ++key_iter, ++value_iter;
    }
}

//---------------------------------------------------------------------------//
// TRIM FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Trim a sorted container by cutting off everything that lies outside
 *        the given bounds.
 */
template<class Container>
void trim(Container                            &container,
          const typename Container::value_type &lower_bound,
          const typename Container::value_type &upper_bound)
{
    typedef typename Container::value_type               val_type;
    typedef typename Compare_Traits<val_type>::eql_comp  eql_cmp;
    typedef typename Compare_Traits<val_type>::less_comp less_cmp;

    REQUIRE( profugus::is_sorted(container.begin(), container.end()) );

    // Trim
    return trim(container, lower_bound, upper_bound, eql_cmp(), less_cmp());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Trim a sorted container by cutting off everything that lies outside
 *        the given bounds within epsilon.
 */
template<class Container>
void trim(Container                            &container,
          const typename Container::value_type &lower_bound,
          const typename Container::value_type &upper_bound,
          const typename Container::value_type &eps)
{
    typedef typename Container::value_type                 value_type;

    REQUIRE( profugus::is_sorted(container.begin(), container.end(),
                                  SoftIsLess<value_type>(eps)) );

    // Trim
    return trim(container, lower_bound, upper_bound,
                SoftIsEqual<value_type>(eps), SoftIsLess<value_type>(eps));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Trim a sorted container so that all of its values lie between
 *        \a lower_bound and \a upper_bound within epsilon
 */
template<class Container, class EqualOp, class LessOp>
void trim(Container                            &container,
          const typename Container::value_type &lower_bound,
          const typename Container::value_type &upper_bound,
          EqualOp                               equal_comp,
          LessOp                                less_comp)
{
    typedef typename Container::iterator      iterator;

    REQUIRE( profugus::is_sorted(container.begin(), container.end(),
                                  less_comp) );
    REQUIRE( less_comp(lower_bound, upper_bound) );

    // Get iterators to the beginning and end of the range we want to *keep*
    iterator lower_iter =
        std::lower_bound(container.begin(), container.end(), lower_bound,
                         less_comp);
    iterator upper_iter =
        std::lower_bound(container.begin(), container.end(), upper_bound,
                         less_comp);

    // Adjust upper iter if it equals the upper bound
    if(upper_iter != container.end() && equal_comp(*upper_iter, upper_bound))
    {
        ++upper_iter;
    }

    // Calculate the distances
    int dist_1 = std::distance(container.begin(), lower_iter);

    // Erase the trimmed end
    container.erase(upper_iter, container.end());
    container.erase(container.begin(), container.begin()+dist_1);
}

//---------------------------------------------------------------------------//
// MAKE_UNIQUE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Make the elements in a container unique.  Does not assume the
 *        container is sorted.
 */
template<class Container>
void make_unique(Container &container)
{
    typedef typename Container::value_type               val_type;
    typedef typename Compare_Traits<val_type>::less_comp less_cmp;
    typedef typename Compare_Traits<val_type>::eql_comp  eql_cmp;

    make_unique(container, eql_cmp(), less_cmp());
}

//---------------------------------------------------------------------------//
// MAKE UNIQUE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Make the elements in a container unique within epsilon.  Does not
 *        assume the container is sorted.
 */
template<class Container>
void make_unique(Container                            &container,
                 const typename Container::value_type &eps)
{
    typedef typename Container::value_type                 value_type;

    make_unique(container, SoftIsEqual<value_type>(eps),
                SoftIsLess<value_type>(eps));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns a container that has been made unique.
 */
template<class Container, class EqualOp, class LessOp>
void make_unique(Container &container,
                 EqualOp    equal_comp,
                 LessOp     less_comp)
{
    // Do a sort
    std::sort(container.begin(), container.end(), less_comp);

    // Make the container unique
    typename Container::iterator iter =
        std::unique(container.begin(), container.end(), equal_comp);

    // Erase the non-unique entries at the end
    container.erase(iter, container.end());
    ENSURE( is_unique(container.begin(), container.end(),
                       equal_comp, less_comp) );
}

//---------------------------------------------------------------------------//
// NUMERICAL FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return the arithmetic mean of the elements in a container.
 */
template<class InputIterator>
typename std::iterator_traits<InputIterator>::value_type
mean(InputIterator begin,
     InputIterator end)
{
    typedef typename ITT<InputIterator>::vt  value_type;

    // Accumulate the container elements
    value_type ret = std::accumulate(begin, end, value_type(0));
    return ret / std::distance(begin, end);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the median of the elements in a container.
 */
template<class InputIterator>
typename std::iterator_traits<InputIterator>::value_type
median(InputIterator begin,
       InputIterator end)
{
    typedef typename ITT<InputIterator>::vt    value_type;

    // Create a vector to hold the information
    std::vector<value_type> vec( std::distance(begin, end) );

    // Copy and sort the elements
    std::partial_sort_copy(begin, end, vec.begin(), vec.end());

    // Find the median element
    if(vec.size() % 2 == 0)
    {
        // Even number of elements.  The middle element is the average between
        // the two
        value_type val_1 = vec[vec.size()/2 - 1];
        value_type val_2 = vec[vec.size()/2];
        return (val_1 + val_2) / 2;
    }
    else
    {
        // Odd number of elements.  Return the middle element
        return vec[vec.size()/2];
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Truncate a vector to [beg, end), renormalizing it to 1
 *
 * We don't allow beg = end because a zero-length vector can't be normalized.
 * For obvious reasons, this should only be used for floating point values.
 *
 * \return Normalization constant: sum of original entries between [beg, end)
 */
template<class Container>
typename Container::value_type
truncate_and_norm(typename Container::iterator  first,
                  typename Container::iterator  last,
                  Container                    &cont)
{
    typedef typename Container::iterator       iterator;
    typedef typename Container::const_iterator const_iterator;
    typedef typename Container::value_type     value_type;

    REQUIRE(std::distance(cont.begin(), first) >= 0);
    REQUIRE(std::distance(last, cont.end()) >= 0);
    REQUIRE(std::distance(first, last) > 0);

    // Calculate the sum of the values in the truncated range
    const value_type norm = std::accumulate(first, last, value_type(0));

    // Copy normalized range. We can do this inline because src >= dst
    iterator dst = cont.begin();
    for (const_iterator src = first; src != last; ++src, ++dst)
    {
        *dst = *src / norm;
    }

    // Delete the remainder of the container
    cont.erase(dst, cont.end());

    ENSURE(cont.size() == static_cast<unsigned int>(last - first));
    ENSURE(profugus::soft_equiv(
                value_type(1),
                std::accumulate(cont.begin(), cont.end(), value_type(0))));
    return norm;
}

} // end namespace profugus

#endif // utils_Container_Functions_i_hh

//---------------------------------------------------------------------------//
//              end of utils/Container_Functions.i.hh
//---------------------------------------------------------------------------//
