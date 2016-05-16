//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/Flat_Table.hh
 * \author Seth R Johnson
 * \date   Fri Apr 11 15:42:54 2014
 * \brief  Flat_Table class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Flat_Table_hh
#define Utils_utils_Flat_Table_hh

#include <utility>
#include <vector>

#include "View_Field.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Flat_Table
 * \brief Linearized contiguous table with ragged right edges
 *
 * The flat table acts like a "vector of vectors" but has a construction phase
 * and completion phase. It is constructed and accessed in row-major format: a
 * single row is stored contiguously. This provides better cache/memory
 * coherency when searching over the table, and it reduces the storage
 * requirements by requiring (in the limit of many rows) only a single
 * size_type per row rather than three pointers per row.
 *
 * Construction occurs by creating new rows (thus preventing the prior row's
 * size from changing) and appending to the current row. Two \c reserve
 * methods will help obviate the performance/memory penalty of repeated \c
 * push_back and \c insert operations, and a \c shrink_to_fit method will
 * reallocate the table using just enough space to hold the elements:
 * \code
    typedef Flat_Table<double>             Flat_Table_t;
    typedef Flat_Table_t::const_View_Field const_View_Field;
    Flat_Table_t ft;

    // First row: [4.125, 3]
    ft.start_row();
    ft.extend_row(4.125);
    ft.extend_row(3);
    // Second row: [2]
    ft.start_row();
    ft.extend_row(2);
    // Third row: [] (empty)
    ft.start_row();
    // Fourth row: [1,2,3,4]
    int vals[] = {1,2,3,4};

    ft.start_row();
    ft.insert(std::begin(vals), std::end(vals));
 *
 * To build a flat table:
 * \code
    typedef Flat_Table<double> Table_t;
    Table_t tab;
    tab.start_row();
    tab.extend_row(1.23);
    tab.start_row();
    tab.extend_row(some_vec.begin(), some_vec.end());
   \endcode
 *
 * Because we don't have an iterator class that acts like a pointer to
 * View_Field, we don't currently support iteration over this container.
 * Instead, loop over indices and extract views:
 * \code
    for (std::size_t j = 0; j < ft.size(); ++j)
    {
        const_View_Field vf(ft[i]);
        for (std::size_t i = 0; i < vf.size(); ++i)
        {
            vf[i] = i + j;
        }
    }
   \endcode
 *
 * Note that adding to the table in most circumstances will invalidate all
 * iterators and views to the data.
 */
/*!
 * \example utils/test/tstFlat_Table.cc
 *
 * Test of Flat_Table.
 */
//===========================================================================//

template <class T>
class Flat_Table
{
    typedef Flat_Table<T> This;
  public:
    typedef T element_type;

    //@{
    //! Container typedefs.
    typedef std::vector<element_type>             Container_T;
    typedef typename Container_T::size_type       size_type;
    typedef typename Container_T::reference       reference;
    typedef typename Container_T::const_reference const_reference;
    //@}

    typedef const_View_Field<element_type> const_View_Field_t;
    typedef View_Field<element_type>       View_Field_t;

    typedef const_View_Field<size_type> const_View_Field_Index;

    typedef std::pair<const_reference, const_reference> pair_cref;

  private:
    // >>> DATA
    typedef std::vector<size_type> Vec_Index;

    //! Linearized items
    Container_T d_storage;
    //! Begin index (into d_storage) of row i of data
    Vec_Index d_row_index;

  public:

    // >>> CONSTRUCTION

    // Create a new table
    Flat_Table();

    // Create a new table explicitly from offsets/values
    Flat_Table(const_View_Field_t values, const_View_Field_Index offsets);

    // Create a new table from a 2D field
    template<class InputIteratorIterator>
    inline Flat_Table(InputIteratorIterator first, InputIteratorIterator last);

    //! Create a new table from an initializer list
    Flat_Table(std::initializer_list<std::initializer_list<T>> init)
        : Flat_Table(init.begin(), init.end())
    { /* * */ }

    // Reserve the number of sub-elements (not the number of rows) expected
    inline void reserve(size_type total_size);

    // Reserve the number of rows expected
    inline void reserve_rows(size_type size);

    // Begin a new row
    inline void start_row();

    // Extend the current row
    template<class InputIterator>
    inline void extend_row(InputIterator first, InputIterator last);

    // Extend the current row with a single element
    inline void extend_row(const_reference element);

    // Shrink the storage and index vectors by copying to a new table
    inline void shrink_to_fit();

    // Remove all rows, all columns
    inline void clear();

    // Swap with another one of us (cheap)
    inline void swap(This& rhs);

    // Append another table to the end of ours
    inline void extend(const This& rhs);

    // >>> ACCESSORS

    //! Number of rows
    size_type num_rows() const { return d_row_index.size() - 1; }

    //! Number of table elements
    size_type num_elements() const { return d_storage.size(); }

    //! Whether no rows have been added (no values either)
    size_type empty() const { return num_rows() == 0; }

    //! Number of rows
    size_type size() const { return num_rows(); }

    //! Number of total elements stored
    size_type total_size() const { return num_elements(); }

    // Whether the table is rectangular (non-jagged)
    inline bool is_rectangular() const;

    // >>> ROW ACCESSORS

    // Access a row of data
    inline const_View_Field_t row(size_type row_index) const;

    //! Access a view field of a row by row index
    const_View_Field_t operator[](size_type row_index) const
    {
        REQUIRE(row_index < size());
        return row(row_index);
    }

    // Access a row of data
    inline View_Field_t row(size_type row_index);

    //! Access a view field of a row by row index
    View_Field_t operator[](size_type row_index)
    {
        REQUIRE(row_index < size());
        return row(row_index);
    }

    // Access a single element
    inline reference get(size_type row_index, size_type col_index);

    // Access a single element
    inline const_reference get(
            size_type row_index,
            size_type col_index) const;

    // >>> ADVANCED FUNCTIONALITY

    // Get the element, and the next one in the storage vector.
    inline pair_cref get_with_next(
            size_type row_index,
            size_type col_index) const;

    //! Get a view into the flattened data.
    const_View_Field_t values() const
    {
        return profugus::make_view(d_storage);
    }

    //! Get a view into the storage offset of each row.
    const_View_Field_Index offsets() const
    {
        return profugus::make_view(d_row_index);
    }

  protected:
    // Get the index into d_storage from a row and column.
    inline size_type index_(size_type row, size_type col) const;

};

//---------------------------------------------------------------------------//
// INLINE TEMPLATE FUNCTIONS
//---------------------------------------------------------------------------//
//! Cheap swap operation
template <class T>
inline void swap(Flat_Table<T>& a, Flat_Table<T>& b) { a.swap(b); }

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE TEMPLATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

#include "Flat_Table.i.hh"

//---------------------------------------------------------------------------//
#endif // Utils_utils_Flat_Table_hh

//---------------------------------------------------------------------------//
//                 end of Flat_Table.hh
//---------------------------------------------------------------------------//
