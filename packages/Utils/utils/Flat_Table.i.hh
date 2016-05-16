//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/Flat_Table.i.hh
 * \author Seth R Johnson
 * \date   Fri Apr 11 15:42:54 2014
 * \brief  Member definitions of class Flat_Table.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Flat_Table_i_hh
#define Utils_utils_Flat_Table_i_hh

#include <iterator>
#include <algorithm>

namespace profugus
{
//---------------------------------------------------------------------------//
// CONSTRUCTION
//---------------------------------------------------------------------------//
/*!
 * \brief Create an empty table
 *
 * This initializes the table with no rows and no columns, analogous to an
 * empty vector of vectors.
 */
template<class T>
Flat_Table<T>::Flat_Table()
{
    clear();

    ENSURE(d_row_index.size() == 1);
    ENSURE(size() == 0);
    ENSURE(total_size() == 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a new table explicitly from values and offsets.
 *
 * This is primarily used for reading tables from HDF5.
 */
template<class T>
Flat_Table<T>::Flat_Table(
        const_View_Field_t      values,
        const_View_Field_Index  offsets)
    : d_storage(values.begin(), values.end())
    , d_row_index(offsets.begin(), offsets.end())
{
    REQUIRE(!offsets.empty());
    REQUIRE(offsets.front() == 0u);
    REQUIRE(offsets.back() == values.size());
    REQUIRE(std::is_sorted(offsets.begin(), offsets.end()));

    ENSURE(!d_row_index.empty());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a new table from a 2D field.
 */
template<class T>
template<class InputIteratorIterator>
Flat_Table<T>::Flat_Table(
        InputIteratorIterator row_iter,
        InputIteratorIterator row_last)
{
    clear();

    for (; row_iter != row_last; ++row_iter)
    {
        start_row();
        extend_row(row_iter->begin(), row_iter->end());
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reserve the number of sub-elements (not the number of rows)
 */
template<class T>
void Flat_Table<T>::reserve(size_type total_size)
{
    d_storage.reserve(total_size);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reserve the number of rows expected
 */
template<class T>
void Flat_Table<T>::reserve_rows(size_type size)
{
    d_row_index.reserve(size + 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a new row, ending the previous row
 */
template<class T>
void Flat_Table<T>::start_row()
{
    REMEMBER(size_type orig_num_rows = num_rows());

    d_row_index.push_back(num_elements());

    ENSURE(num_rows() == orig_num_rows + 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Extend the current row by the given elements
 */
template<class T>
template<class InputIterator>
void Flat_Table<T>::extend_row(InputIterator first, InputIterator last)
{
    REQUIRE(num_rows() > 0); // start_row must be called first
    CHECK(!d_row_index.empty());

    d_storage.insert(d_storage.end(), first, last);
    d_row_index.back() = num_elements();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Extend the current row with a single element
 */
template<class T>
void Flat_Table<T>::extend_row(const_reference element)
{
    REQUIRE(num_rows() > 0); // start_row must be called first
    CHECK(!d_row_index.empty());

    d_storage.push_back(element);
    d_row_index.back() = d_storage.size();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Shrink the storage and index vectors by copying to a new table
 */
template<class T>
void Flat_Table<T>::shrink_to_fit()
{
    using std::swap;
    // Allocate new properly sized data vector and swap
    {
        Container_T new_data(d_storage);
        swap(new_data, d_storage);
    }

    // Allocate new properly sized index vector and swap
    {
        Vec_Index new_indices(d_row_index);
        swap(new_indices, d_row_index);
    }

    ENSURE(d_storage.capacity() == d_storage.size());
    ENSURE(d_row_index.capacity() == d_row_index.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reset: zero rows with zero elements
 */
template<class T>
void Flat_Table<T>::clear()
{
    d_row_index.assign(1,0);
    d_storage.clear();
    ENSURE(size() == 0);
    ENSURE(total_size() == 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Swap with another one of us (cheap)
 */
template<class T>
void Flat_Table<T>::swap(This& rhs)
{
    if (this == &rhs)
        return;

    d_row_index.swap(rhs.d_row_index);
    d_storage.swap(rhs.d_storage);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Extend with the contents of another table
 */
template<class T>
void Flat_Table<T>::extend(const This& other)
{
    REQUIRE(this != &other);

    const size_type original_size = d_storage.size();
    CHECK(original_size == d_row_index.back());

    // Append elements
    d_storage.insert(d_storage.end(),
                     other.d_storage.begin(), other.d_storage.end());

    // Append rows
    d_row_index.reserve(d_row_index.size() + other.num_rows());
    auto idx_iter = other.d_row_index.begin() + 1;
    for (auto end_iter = other.d_row_index.end(); idx_iter != end_iter;
         ++idx_iter)
    {
        d_row_index.push_back(original_size + *idx_iter);
    }
}

//---------------------------------------------------------------------------//
// ACCESSORS
//---------------------------------------------------------------------------//
/*!
 * \brief Access view for a single row
 */
template<class T>
typename Flat_Table<T>::const_View_Field_t
Flat_Table<T>::row(size_type row_index) const
{
    REQUIRE(row_index < size());
    const size_type start = d_row_index[row_index];
    const size_type stop  = d_row_index[row_index + 1];

    ENSURE(start <= stop);
    ENSURE(stop <= d_storage.size());
    return const_View_Field_t(d_storage.data() + start,
                              d_storage.data() + stop);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access view for a single row
 */
template<class T>
typename Flat_Table<T>::View_Field_t
Flat_Table<T>::row(size_type row_index)
{
    REQUIRE(row_index < size());
    const size_type start = d_row_index[row_index];
    const size_type stop  = d_row_index[row_index + 1];

    ENSURE(start <= stop);
    ENSURE(stop <= d_storage.size());
    return View_Field_t(d_storage.data() + start, d_storage.data() + stop);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a single element
 */
template<class T>
typename Flat_Table<T>::const_reference
Flat_Table<T>::get(size_type row_index, size_type col_index) const
{
    return d_storage[index_(row_index, col_index)];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a single element
 */
template<class T>
typename Flat_Table<T>::reference
Flat_Table<T>::get(size_type row_index, size_type col_index)
{
    return d_storage[index_(row_index, col_index)];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the element, and the next one in the storage vector.
 *
 * This is used for creating panels. You may not attempt to access the final
 * stored entry.
 */
template<class T>
typename Flat_Table<T>::pair_cref
Flat_Table<T>::get_with_next(size_type row_index, size_type col_index) const
{
    const size_type first_index = index_(row_index, col_index);
    REQUIRE(first_index + 1 < d_storage.size());

    return pair_cref(d_storage[first_index], d_storage[first_index + 1]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Whether the table is rectangular (all rows have same size).
 *
 * If this returns true, it would be better to have a single vector with a
 * Hyperslab view into it: that structure requires one fewer memory
 * dereference.
 */
template<class T>
bool Flat_Table<T>::is_rectangular() const
{
    auto iter = d_row_index.begin();
    CHECK(iter != d_row_index.end());

    // Get the first row size
    size_type prev_index = *iter;
    ++iter;
    const size_type first_row_size = *iter - prev_index;

    for (auto end_iter = d_row_index.end(); iter != end_iter; ++iter)
    {
        const size_type cur_index = *iter;
        if ((cur_index - prev_index) != first_row_size)
            return false;
        prev_index = cur_index;
    }

    return true;
}

//---------------------------------------------------------------------------//
// PROTECTED
//---------------------------------------------------------------------------//
template<class T>
typename Flat_Table<T>::size_type
Flat_Table<T>::index_(size_type row, size_type col) const
{
    REQUIRE(row < num_rows());
    const size_type begin = d_row_index[row];
    REQUIRE(begin + col < d_row_index[row + 1]);

    ENSURE(begin + col < d_storage.size());
    return begin + col;
}

} // end namespace profugus

#endif // Utils_utils_Flat_Table_i_hh

//---------------------------------------------------------------------------//
//                 end of Flat_Table.i.hh
//---------------------------------------------------------------------------//
