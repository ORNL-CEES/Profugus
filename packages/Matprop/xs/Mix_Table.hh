//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Matprop/xs/Mix_Table.hh
 * \author Seth R Johnson
 * \date   Thu Feb 04 15:02:21 2016
 * \brief  Mix_Table class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Matprop_xs_Mix_Table_hh
#define Matprop_xs_Mix_Table_hh

#include "Utils/utils/Flat_Table.hh"
#include "Utils/utils/Flat_Table.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Mix_Table
 * \brief Sparse mix table: rows are new matids, cols are unmixed matids
 *
 * The mix table is a sparse matrix whose rows are volume fractions of
 * original, "unmixed" materials.
 */
/*!
 * \example material/test/tstMix_Table.cc
 *
 * Test of Mix_Table.
 */
//===========================================================================//

class Mix_Table
{
  public:
    //@{
    //! Typedefs
    typedef double                                 value_type;
    typedef int                                    matid_type;
    typedef std::pair<matid_type, value_type>      Matid_Frac;
    typedef profugus::const_View_Field<Matid_Frac>  const_View_Matfrac;
    typedef profugus::const_View_Field<matid_type>  const_View_Mat;
    typedef profugus::const_View_Field<value_type>  const_View_Frac;
    typedef profugus::const_View_Field<std::size_t> const_View_Idx;
    typedef profugus::View_Field<matid_type>        View_Mat;
    typedef profugus::View_Field<value_type>        View_Frac;
    //@}

  private:
    // >>> DATA

    //! Table storage
    profugus::Flat_Table<Matid_Frac> d_table;

    //! Highest unmixed ID
    matid_type d_num_cols;

  public:

    // Build from sparse coordinates
    static Mix_Table from_coo(const_View_Mat  rows,
                              const_View_Mat  cols,
                              const_View_Frac values);

  public:

    // >>> CONSTRUCTION

    // Empty constructor
    Mix_Table();

#ifndef SWIG
    // Default copy constructor
    Mix_Table(const Mix_Table&) = default;
    Mix_Table& operator=(const Mix_Table&) = default;

    // Un-complete other table on move
    inline Mix_Table(Mix_Table&&);
    inline Mix_Table& operator=(Mix_Table&&);
#endif

    // Begin a new row (mixed matid)
    void start_row() { d_table.start_row(); }

    // Add a pure matid to this row
    inline void extend_row(matid_type unmixed_id, value_type fraction);

    // Remove all rows, all columns
    inline void clear();

    // Normalize all rows, update maximum column
    void complete();

    // Swap with another one of us (cheap)
    inline void swap(Mix_Table& rhs);

    // >>> ACCESSORS

    //! Whether the mix table has been completed
    bool completed() const { return d_num_cols != 0; }

    //! Number of rows (1 + highest mixed matid)
    matid_type num_rows() const { return d_table.num_rows(); }

    //! Number of cols (1 + highest unmixed matid)
    matid_type num_cols() const { REQUIRE(completed()); return d_num_cols; }

    //! Number of nonzero elements
    matid_type num_elements() const { return d_table.num_elements(); }

    //! Whether no rows have been added (no values either)
    bool empty() const { return d_table.empty(); }

    // For a given row, access nonzero column/fraction values
    inline const_View_Matfrac row(matid_type row_id) const;

    //! For a given row, access nonzero column/fraction values
    const_View_Matfrac operator[](matid_type mixed_id) const
    {
        return row(mixed_id);
    }

    // Access a particular row/column value
    value_type find_fraction(matid_type mixed_id, matid_type pure_id) const;

    // Concatenate another mix table onto the end of ours
    void extend(const Mix_Table& additional);

    // >>> FLATTENED ACCESS

    // Get the COO representation of this table
    void to_coo(View_Mat row, View_Mat col, View_Frac val) const;

    // Return a flattened list of all values
    const_View_Frac values() const;

    // Return a flattened list of all columns
    const_View_Mat columns() const;

    // Return a flattened list of all row start indices in the above lists
    const_View_Idx row_offsets() const;
};

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//
//! Cheap swap operation
inline void swap(Mix_Table& a, Mix_Table& b) { a.swap(b); }

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Mix_Table.i.hh"

//---------------------------------------------------------------------------//
#endif // Matprop_xs_Mix_Table_hh

//---------------------------------------------------------------------------//
// end of Transcore/material/Mix_Table.hh
//---------------------------------------------------------------------------//
