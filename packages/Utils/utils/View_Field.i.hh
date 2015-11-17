//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/View_Field.i.hh
 * \author Gregory G. Davidson
 * \date   Tue Feb 16 08:52:06 2010
 * \brief  Member definitions of View_Field and const_View_Field classes.
 * \note   Copyright (C) 2010 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_View_Field_i_hh
#define Utils_utils_View_Field_i_hh

#include <cstring>

namespace profugus
{

//===========================================================================//
// View_Field INLINE DEFINITIONS
//===========================================================================//
/*!
 * \brief Default constructor.  Creates an empty view.
 */
template<typename T>
View_Field<T>::View_Field()
{
    /* empty */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor. Takes two pointers and an optional stride parameter.
 */
template<typename T>
View_Field<T>::View_Field(pointer     first,
                          pointer     last,
                          stride_type stride)
    : d_begin_iterator(first, stride)
    , d_end_iterator(last, stride)
{
    REQUIRE(std::distance(first, last) >= 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct from another set of view field iterators
 */
template<typename T>
View_Field<T>::View_Field(iterator    first,
                          iterator    last)
    : d_begin_iterator(first)
    , d_end_iterator(last)
{
    REQUIRE(std::distance(first, last) >= 0);
    REQUIRE(first.stride() == last.stride());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign function from another View_Field
 *
 * This function is optimized for stride 1 Fiew_Fields.
 */
template<typename T>
void View_Field<T>::assign(const_iterator first,
                           const_iterator last)
{
    REQUIRE(last - first == d_end_iterator - d_begin_iterator);
    REQUIRE(first.stride() == last.stride());

    // Optimization for stride == 1
    if (first.stride() == 1 && d_begin_iterator.stride() == 1)
    {
        // Static assert not compiling to trivial types...?
        std::memmove(d_begin_iterator.get_pointer(), first.get_pointer(),
                     sizeof(T) * (d_end_iterator-d_begin_iterator));
    }
    else
    {
        std::copy(first, last, d_begin_iterator);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief General-purpose assign function.
 */
template<typename T>
template<typename InputIterator>
void View_Field<T>::assign(InputIterator first,
                           InputIterator last)
{
    REQUIRE(last - first == size());

    std::copy(first, last, d_begin_iterator);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns data reference at index \a i.
 */
template<typename T>
typename View_Field<T>::reference
View_Field<T>::operator[](size_type i)
{
    REQUIRE(valid_index(i));

    return *(d_begin_iterator + i);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns const data reference at index \a i.
 */
template<typename T>
typename View_Field<T>::const_reference
View_Field<T>::operator[](size_type i) const
{
    REQUIRE(valid_index(i));

    return *(d_begin_iterator + i);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the front piece of data.
 */
template<typename T>
typename View_Field<T>::reference
View_Field<T>::front()
{
    REQUIRE(!empty());

    return *d_begin_iterator;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the front piece of data.
 */
template<typename T>
typename View_Field<T>::const_reference
View_Field<T>::front() const
{
    REQUIRE(!empty());

    return *d_begin_iterator;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the back piece of data.
 */
template<typename T>
typename View_Field<T>::reference
View_Field<T>::back()
{
    REQUIRE(!empty());

    return *(d_end_iterator-1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the back piece of data.
 */
template<typename T>
typename View_Field<T>::const_reference
View_Field<T>::back() const
{
    REQUIRE(!empty());

    return *(d_end_iterator-1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the number of field values.
 */
template<typename T>
typename View_Field<T>::size_type
View_Field<T>::size() const
{
    return std::distance(d_begin_iterator, d_end_iterator);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the stride between successive elements.
 */
template<typename T>
typename View_Field<T>::stride_type
View_Field<T>::stride() const
{
    ENSURE(d_begin_iterator.stride() == d_end_iterator.stride());
    return d_begin_iterator.stride();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Obtain raw access to underlying data.
 *
 * \note That will throw an error when the view field is strided! This is to
 * prevent unintentionally trying to send a strided data field to a
 * command that assumes contiguous data when it gets a pointer.
 */
template<typename T>
typename View_Field<T>::pointer
View_Field<T>::data()
{
    REQUIRE(stride() == 1);
    return d_begin_iterator.get_pointer();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Obtain raw access to underlying data.
 */
template<typename T>
typename View_Field<T>::const_pointer
View_Field<T>::data() const
{
    REQUIRE(stride() == 1);
    return d_begin_iterator.get_pointer();
}

//---------------------------------------------------------------------------//
/*!
 * \brief CHECKs that index \a i is within bounds.
 */
template<typename T>
bool View_Field<T>::valid_index(size_type i) const
{
    return d_begin_iterator + i < d_end_iterator;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Begin reverse const_iterator.
 */
template<typename T>
typename View_Field<T>::const_reverse_iterator
View_Field<T>::rbegin() const
{
    return const_reverse_iterator(end());
}

//---------------------------------------------------------------------------//
/*!
 * \brief End reverse const_iterator.
 */
template<typename T>
typename View_Field<T>::const_reverse_iterator
View_Field<T>::rend() const
{
    return const_reverse_iterator(begin());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Divides the field into equal slices of size \a slice_size, and
 *        returns the \a slice_num-th slice (starting from 0).
 */
template<typename T>
View_Field<T> View_Field<T>::slice(size_type slice_size,
                                   int       slice_num) const
{
    REQUIRE(!empty());
    REQUIRE(stride() == 1);

    REQUIRE(slice_size > 0);
    REQUIRE(slice_num >= 0);
    REQUIRE(slice_size <= size());
    REQUIRE(slice_num * slice_size <= size());

    // Calculate the begin and end indexes
    size_type index_begin = slice_num * slice_size;
    size_type index_end   = index_begin + slice_size;
    CHECK(d_begin_iterator + index_end <= d_end_iterator);

    // Create and return a View_Field
    return View_Field_t(d_begin_iterator.get_pointer() + index_begin,
                        d_begin_iterator.get_pointer() + index_end);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Divides the field into equal slices of size \a slice_size, and
 *        returns slices slice_begin to slice_end.
 */
template<typename T>
View_Field<T> View_Field<T>::slice(size_type slice_size,
                                   int       slice_begin,
                                   int       slice_end) const
{
    REQUIRE(!empty());
    REQUIRE(stride() == 1);

    REQUIRE(slice_size > 0);
    REQUIRE(slice_begin >= 0);
    REQUIRE(slice_end >= slice_begin);
    REQUIRE(slice_size * slice_begin <= size());
    REQUIRE(slice_size * slice_end <= size());

    // Calculate the begin and end indexes
    size_type index_begin = slice_size * slice_begin;
    size_type index_end   = slice_size * slice_end;
    CHECK(d_begin_iterator + index_end <= d_end_iterator);

    // Create and return a View_Field
    return View_Field_t(d_begin_iterator.get_pointer() + index_begin,
                        d_begin_iterator.get_pointer() + index_end);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns a slice of the field from (\a index_begin, \a index_end].
 */
template<typename T>
View_Field<T> View_Field<T>::general_slice(size_type index_begin,
                                           size_type index_end) const
{
    REQUIRE(d_begin_iterator);
    REQUIRE(d_end_iterator);
    REQUIRE(!empty());

    REQUIRE(index_end >= index_begin);
    REQUIRE(d_begin_iterator + index_begin < d_end_iterator);
    REQUIRE(d_begin_iterator + index_end <= d_end_iterator);

    return View_Field_t(
            d_begin_iterator.get_pointer() + index_begin,
            d_begin_iterator.get_pointer() + index_end * stride(),
            stride());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns a strided slice of the view field.
 */
template<typename T>
View_Field<T> View_Field<T>::strided_slice(
        size_type start,
        size_type stop,
        size_type step) const
{
    REQUIRE(start <= stop);
    REQUIRE(step > 0);
    REQUIRE(stop <= size() + step - 1);
    REQUIRE((stop - start) % step == 0);

    return View_Field_t(d_begin_iterator.get_pointer() + start*stride(),
                        d_begin_iterator.get_pointer() + stop*stride(),
                        stride()*step);
}

//===========================================================================//
// CONST_FIELD_VIEW INLINE DEFINITIONS
//===========================================================================//
/*!
 * \brief Default constructor.  Creates an empty view.
 */
template<typename T>
const_View_Field<T>::const_View_Field()
{
    /* empty */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor. Build from a View_Field.
 */
template<typename T>
const_View_Field<T>::const_View_Field(const View_Field_t& rhs)
    : d_begin_iterator(rhs.begin())
    , d_end_iterator(rhs.end())
{
    ENSURE(size() == rhs.size());
    ENSURE(stride() == rhs.stride());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor. Takes two pointers and an optional stride.
 */
template<typename T>
const_View_Field<T>::const_View_Field(const_pointer first,
                                      const_pointer last,
                                      stride_type   stride)
    : d_begin_iterator(first, stride)
    , d_end_iterator(last, stride)
{
    REQUIRE(std::distance(first, last) >= 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct from another set of view field iterators
 */
template<typename T>
const_View_Field<T>::const_View_Field(const_iterator    first,
                                      const_iterator    last)
    : d_begin_iterator(first)
    , d_end_iterator(last)
{
    REQUIRE(std::distance(first, last) >= 0);
    REQUIRE(first.stride() == last.stride());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns constant field value at index i.
 */
template<typename T>
typename const_View_Field<T>::const_reference
const_View_Field<T>::operator[](size_type i) const
{
    REQUIRE(valid_index(i));

    return *(d_begin_iterator + i);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the piece of data at the front of the container.
 */
template<typename T>
typename const_View_Field<T>::const_reference
const_View_Field<T>::front() const
{
    REQUIRE(!empty());

    return *d_begin_iterator;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the piece of data at the back of the container.
 */
template<typename T>
typename const_View_Field<T>::const_reference
const_View_Field<T>::back() const
{
    REQUIRE(!empty());

    return *(d_end_iterator - 1);
}


//---------------------------------------------------------------------------//
/*!
 * \brief Return the number of field values.
 */
template<typename T>
typename const_View_Field<T>::size_type
const_View_Field<T>::size() const
{
    return std::distance(d_begin_iterator, d_end_iterator);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the stride between successive elements.
 */
template<typename T>
typename const_View_Field<T>::stride_type
const_View_Field<T>::stride() const
{
    ENSURE(d_begin_iterator.stride() == d_end_iterator.stride());
    return d_begin_iterator.stride();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Obtain raw access to underlying data.
 */
template<typename T>
typename const_View_Field<T>::const_pointer
const_View_Field<T>::data() const
{
    REQUIRE(stride() == 1);
    return d_begin_iterator.get_pointer();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Verifies that index \a i is a valid index.
 */
template<typename T>
bool const_View_Field<T>::valid_index(size_type i) const
{
    if (!d_begin_iterator.is_valid())
        return false;
    if (!d_end_iterator.is_valid())
        return false;

    return d_begin_iterator + i < d_end_iterator;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a begin reverse const_iterator.
 */
template<typename T>
typename const_View_Field<T>::const_reverse_iterator
const_View_Field<T>::rbegin() const
{
    return const_reverse_iterator(end());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return an end reverse const_iterator.
 */
template<typename T>
typename const_View_Field<T>::const_reverse_iterator
const_View_Field<T>::rend() const
{
    return const_reverse_iterator(begin());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Divides the field into equal slices of size \a slice_size, and
 *        returns the slice_num-th slice (starting from 0).
 */
template<typename T>
const_View_Field<T> const_View_Field<T>::slice(size_type slice_size,
                                               int       slice_num) const
{
    REQUIRE(!empty());
    REQUIRE(stride() == 1);

    REQUIRE(slice_size > 0);
    REQUIRE(slice_num >= 0);
    REQUIRE(slice_size <= size());
    REQUIRE(slice_num * slice_size <= size());

    // Calculate the begin and end indexes
    size_type index_begin = slice_num * slice_size;
    size_type index_end   = index_begin + slice_size;
    CHECK(d_begin_iterator + index_end <= d_end_iterator);

    // Create and return a View_Field
    return const_View_Field_t(d_begin_iterator.get_pointer() + index_begin,
                              d_begin_iterator.get_pointer() + index_end);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Divides the field into equal slices of size \a slice_size, and
 *        returns slices slice_begin to slice_end.
 */
template<typename T>
const_View_Field<T> const_View_Field<T>::slice(size_type slice_size,
                                               int       slice_begin,
                                               int       slice_end) const
{
    REQUIRE(stride() == 1);
    REQUIRE(slice_size > 0);
    REQUIRE(slice_begin >= 0);
    REQUIRE(slice_end >= slice_begin);
    REQUIRE(slice_size * slice_begin <= size());
    REQUIRE(slice_size * slice_end <= size());

    // Calculate the begin and end indexes
    size_type index_begin = slice_size * slice_begin;
    size_type index_end   = slice_size * slice_end;
    CHECK(d_begin_iterator + index_end <= d_end_iterator);

    // Create and return a View_Field
    return const_View_Field_t(d_begin_iterator.get_pointer() + index_begin,
                              d_begin_iterator.get_pointer() + index_end);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns a slice of the field from (\a index_begin, \a index_end].
 */
template<typename T>
const_View_Field<T>
const_View_Field<T>::general_slice(size_type index_begin,
                                   size_type index_end) const
{
    REQUIRE(index_begin <= index_end);
    REQUIRE(index_end <= size());

    return const_View_Field_t(
            d_begin_iterator.get_pointer() + index_begin,
            d_begin_iterator.get_pointer() + index_end * stride(),
            stride());
}


//---------------------------------------------------------------------------//
/*!
 * \brief Returns a strided slice of the view field.
 */
template<typename T>
const_View_Field<T> const_View_Field<T>::strided_slice(
        size_type start,
        size_type stop,
        size_type step) const
{
    REQUIRE(start <= stop);
    REQUIRE(step > 0);
    REQUIRE(stop <= size() + step - 1);
    REQUIRE((stop - start) % step == 0);

    return const_View_Field_t(d_begin_iterator.get_pointer() + start * stride(),
                              d_begin_iterator.get_pointer() + stop * stride(),
                              stride()*step);
}

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Wraps a C array in a View_Field.
 *
 * \code
 *
 *   double ref[] = {1,2,3};
 *   const_View_Field<double> vf(ref);
 *
 * \endcode
 */
template<typename T, std::size_t N>
View_Field<T> make_view(T (&array)[N])
{
    REQUIRE(N > 0);

    return View_Field<T>(array, array + N);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wraps a C array in a const_View_Field.
 */
template<typename T, std::size_t N>
const_View_Field<T> make_view(const T (&array)[N])
{
    REQUIRE(N > 0);

    return const_View_Field<T>(array, array + N);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wraps a C array in a const_View_Field.
 */
template<typename T, std::size_t N>
const_View_Field<T> make_const_view(const T (&array)[N])
{
    return make_view<T,N>(array);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wraps a std::vector in a const_View_Field.
 */
template<typename T, class A>
View_Field<T> make_view(std::vector<T,A> &vec)
{
    // Empty vector gets an empty view field (data() may be undefined)
    if (vec.empty())
        return View_Field<T>();

    return View_Field<T>(vec.data(), vec.data() + vec.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wraps a std::vector in a const_View_Field.
 */
template<typename T, class A>
const_View_Field<T> make_view(const std::vector<T,A> &vec)
{
    // Empty vector gets an empty view field (data() may be undefined)
    if (vec.empty())
        return const_View_Field<T>();

    return const_View_Field<T>(vec.data(), vec.data() + vec.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wraps a std::vector in a const_View_Field.
 */
template<typename T, class A>
const_View_Field<T> make_const_view(const std::vector<T,A> &vec)
{
    return make_view<T,A>(vec);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Re-emit a view
 */
template<typename T>
View_Field<T> make_view(profugus::View_Field<T> view)
{
    return view;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Re-emit a view
 */
template<typename T>
const_View_Field<T> make_view(profugus::const_View_Field<T> view)
{
    return view;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Re-emit a view
 */
template<typename T>
const_View_Field<T> make_const_view(profugus::const_View_Field<T> view)
{
    return view;
}

}  // end namespace profugus

#endif // Utils_utils_View_Field_i_hh

//---------------------------------------------------------------------------//
//              end of profugus/View_Field.i.hh
//---------------------------------------------------------------------------//
