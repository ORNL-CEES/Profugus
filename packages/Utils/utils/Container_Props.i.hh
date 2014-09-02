//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Container_Props.i.hh
 * \author Gregory G. Davidson
 * \date   Thu Jan 09 11:13:13 2014
 * \brief  Inline function definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Container_Props_i_hh
#define utils_Container_Props_i_hh

#include <vector>

namespace profugus
{

namespace detail
{

// Soft find (needed by is_subset)
template<class InputIterator>
InputIterator soft_find(InputIterator                          begin,
                        InputIterator                          end,
                        const typename ITT<InputIterator>::vt &value,
                        const typename ITT<InputIterator>::vt &eps)
{
    typedef typename ITT<InputIterator>::vt  value_type;

    REQUIRE(std::distance(begin, end) >= 0);

    return std::find_if(begin, end,
                        std::bind2nd(SoftIsEqual<value_type>(eps), value));
}

} // end namespace profugus::detail

//---------------------------------------------------------------------------//
// IS-SORTED FUNCTION
//---------------------------------------------------------------------------//
/*!
 * \brief Return whether a container is sorted.
 */
template<class InputIterator>
bool is_sorted(InputIterator begin,
               InputIterator end)
{
    typedef typename ITT<InputIterator>::vt              val_type;
    typedef typename Compare_Traits<val_type>::less_comp less_cmp;

    return profugus::is_sorted(begin, end, less_cmp());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return whether a container is sorted.
 */
template<class InputIterator, class BinaryOp>
bool is_sorted(InputIterator iter,
               InputIterator last,
               BinaryOp      comp)
{
    REQUIRE(std::distance(iter, last) >= 0);

    if (iter == last)
        return true;

    // Copy the iterator
    InputIterator prev_iter = iter;

    // Increment iter
    ++iter;

    // Loop until we hit the last
    while (iter != last)
    {
        // Compare against the previous value
        if (!comp(*prev_iter, *iter))
            return false;

        // Update
        ++prev_iter;
        ++iter;
    }
    return true;
}

//---------------------------------------------------------------------------//
// IS_SUBSET FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return whether \a cont2 is a subset of \a cont1.
 *
 * \note By subset, we mean that all of the elements in \a cont2 are in \a
 *       cont1.
 * \pre  Both cont1 and cont2 must be sorted.
 */
template<class Container1, class Container2>
bool is_subset(const Container1 &cont1,
               const Container2 &cont2)
{
    // Define the equality operator
    typedef typename Container1::value_type             val_type;
    typedef typename Compare_Traits<val_type>::eql_comp eql_comp;

    REQUIRE( profugus::is_sorted(cont1.begin(), cont1.end()) );
    REQUIRE( profugus::is_sorted(cont2.begin(), cont2.end()) );

    // Create the equal comparator
    eql_comp comp;

    // Keep a search iterator into Container1
    typename Container1::const_iterator iter_1 = cont1.begin();

    // Loop over container2
    for(typename Container2::const_iterator iter_2 = cont2.begin(),
                                        iter_end_2 = cont2.end();
        iter_2 != iter_end_2; ++iter_2)
    {
        // Find *iter_2 in Container 1
        iter_1 = std::find_if(iter_1, cont1.end(), std::bind2nd(comp, *iter_2));

        if(iter_1 == cont1.end()) return false;
    }

    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return whether \a cont2 is a subset of \a cont1 within epsilon
 *        \a epsilon.
 *
 * \note By subset, we mean that all of the elements in \a cont2 are in \a
 *       cont1.
 * \pre  Both cont1 and cont2 must be sorted.
 */
template<class Container1, class Container2>
bool is_subset(const Container1                      &cont1,
               const Container2                      &cont2,
               const typename Container1::value_type &eps)
{
    typedef typename Container1::value_type     value_type;

    REQUIRE( profugus::is_sorted(cont1.begin(), cont1.end(),
                                  SoftIsLessEqual<value_type>(eps)) );
    REQUIRE( profugus::is_sorted(cont2.begin(), cont2.end(),
                                  SoftIsLessEqual<value_type>(eps)) );

    // Create the equal comparator
    SoftIsEqual<value_type> is_eql(eps);

    // Keep a search iterator into Container1
    typename Container1::const_iterator iter_1 = cont1.begin();

    // Loop over container2
    for(typename Container2::const_iterator iter_2 = cont2.begin(),
                                        iter_end_2 = cont2.end();
        iter_2 != iter_end_2; ++iter_2)
    {
        // Find *iter_2 in Container 1
        iter_1 = detail::soft_find(iter_1, cont1.end(), *iter_2, eps);

        if(iter_1 == cont1.end()) return false;
    }

    return true;
}

//---------------------------------------------------------------------------//
// IS_UNIQUE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Returns whether elements in a container are unique.
 */
template<class InputIterator>
bool is_unique(InputIterator begin,
               InputIterator end)
{
    typedef typename ITT<InputIterator>::vt              val_type;
    typedef typename Compare_Traits<val_type>::eql_comp  eql_cmp;
    typedef typename Compare_Traits<val_type>::less_comp less_cmp;

    return is_unique(begin, end, eql_cmp(), less_cmp());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns whether elements in a container are unique within
 *        \a epsilon.
 */
template<class InputIterator>
bool is_unique(InputIterator                          begin,
               InputIterator                          end,
               const typename ITT<InputIterator>::vt &eps)
{
    typedef typename ITT<InputIterator>::vt value_type;

    return is_unique(begin, end,
                     SoftIsEqual<value_type>(eps),
                     SoftIsLess<value_type>(eps));
}


//---------------------------------------------------------------------------//
/*!
 * \brief Returns whether the entries in a container are unique (within
 *        epsilon).
 */
template<class InputIterator, class EqualOp, class LessOp>
bool is_unique(InputIterator begin,
               InputIterator end,
               EqualOp       equal_comp,
               LessOp        less_comp)
{
    typedef typename ITT<InputIterator>::vt value_type;
    typedef std::vector<value_type>         Vec_ValType;

    // Create a destination container
    Vec_ValType cont( std::distance(begin, end) );

    // Sort-and-copy the container into a new container
    std::partial_sort_copy(begin, end, cont.begin(), cont.end(), less_comp);

    // Removed duplicates that are equal within epsilon
    typename Vec_ValType::iterator iter =
        std::unique(cont.begin(), cont.end(), equal_comp);

    return std::distance(begin, end) == std::distance(cont.begin(), iter);
}

//---------------------------------------------------------------------------//
// NUMERICAL FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Check that a range of elements are all positive or zero.
 */
template<class InputIterator>
bool is_positive(InputIterator begin,
                 InputIterator end)
{
    typedef typename ITT<InputIterator>::vt value_type;

    if(std::find_if(begin, end, std::bind2nd(std::less_equal<value_type>(),
                                             value_type(0))) == end)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check that a range of elements are all negative or zero.
 */
template<class InputIterator>
bool is_negative(InputIterator begin,
                 InputIterator end)
{
    typedef typename ITT<InputIterator>::vt value_type;

    if(std::find_if(begin, end, std::bind2nd(std::greater_equal<value_type>(),
                                             value_type(0))) == end)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//---------------------------------------------------------------------------//
/*! \brief Check a container to ensure its contents are non-negative.
 *
 * This will return false if the container is empty.
 */
template<class InputIterator>
bool is_non_negative(InputIterator begin,
                     InputIterator end)
{
    typedef typename ITT<InputIterator>::vt value_type;

    if(std::find_if(begin, end, std::bind2nd(std::less<value_type>(),
                                             value_type(0))) == end)
    {
        return true;
    }
    else
    {
        return false;
    }
}

} // end namespace profugus

#endif // utils_Container_Props_i_hh

//---------------------------------------------------------------------------//
//              end of utils/Container_Props.i.hh
//---------------------------------------------------------------------------//
