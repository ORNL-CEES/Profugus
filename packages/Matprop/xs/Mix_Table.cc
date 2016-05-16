//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Matprop/xs/Mix_Table.cc
 * \author Seth R Johnson
 * \date   Thu Feb 04 15:02:21 2016
 * \brief  Mix_Table class definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Mix_Table.hh"

#include "Utils/utils/Range.hh"
#include "Utils/utils/View_Field_Struct.hh"

namespace profugus
{
//---------------------------------------------------------------------------//
// METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Build from sparse coordinates
 */
Mix_Table Mix_Table::from_coo(const_View_Mat  row_ids,
                              const_View_Mat  col_ids,
                              const_View_Frac values)
{
    VALIDATE(row_ids.size() == values.size(),
             "Row size " << row_ids.size() << " does not match value "
             "size " << values.size());
    VALIDATE(row_ids.size() == col_ids.size(),
             "Row size " << row_ids.size() << " does not match cols "
             "size " << col_ids.size());

    Mix_Table table;

    // Start the first row if not empty
    if (!values.empty())
        table.start_row();

    matid_type row_index = 0;

    auto r = row_ids.begin();
    auto c = col_ids.begin();
    auto v = values.begin();

    for (auto end_r = row_ids.end(); r != end_r; ++r, ++c, ++v)
    {
        while (*r > row_index)
        {
            table.start_row();
            ++row_index;
        }
        table.extend_row(*c, *v);
    }

    table.complete();
    return std::move(table);
}

//---------------------------------------------------------------------------//
// METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Empty constructor
 */
Mix_Table::Mix_Table()
    : d_table()
    , d_num_cols(0)
{
    ENSURE(!completed());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Normalize all rows, update maximum column
 */
void Mix_Table::complete()
{
    REQUIRE(!completed());
    INSIST(num_rows() > 0, "Can't complete an empty mix table");

    matid_type max_col = 0;

    for (auto i : profugus::range(d_table.num_rows()))
    {
        auto row_data = d_table[i];
        if (row_data.empty())
            continue;

        // Calculate normalization constant
        double row_total = 0.0;
        for (Matid_Frac mv : row_data)
        {
            // Get maximum column
            if (mv.first > max_col)
                max_col = mv.first;

            // Update volume fractions
            row_total += mv.second;
        }
        // Should always be positive because "zero" entries are never added
        CHECK(row_total > 0);

        // Normalize row
        const double norm = 1.0 / row_total;
        for (Matid_Frac& mv : row_data)
        {
            mv.second *= norm;
        }

        // Sort row by pure IDs
        std::sort(row_data.begin(), row_data.end(),
                  [](const Matid_Frac& a, const Matid_Frac& b)
                  { return a.first < b.first; });
    }

    d_num_cols = max_col + 1;

    ENSURE(completed());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a particular row/column value.
 *
 * If called before "completed", either row or column can be outside the mix
 * table, and we return zero rather than raising an
 * error. Additionally the fractions will be unnormalized.
 *
 * If called after completed, we require both mixed and pure to be within the
 * bounds of the mix table.
 */
Mix_Table::value_type Mix_Table::find_fraction(
        matid_type mixed_id,
        matid_type pure_id) const
{
    REQUIRE(mixed_id >= 0);
    REQUIRE(pure_id >= 0);

    if (completed())
    {
        REQUIRE(mixed_id < num_rows());
        REQUIRE(pure_id  < num_cols());

        // Find sorted pure ID in the row
        auto row_view = row(mixed_id);
        auto iter = std::lower_bound(row_view.begin(), row_view.end(), pure_id,
                                     [](const Matid_Frac &mf, matid_type m)
                                     { return mf.first < m; });

        if (iter != row_view.end() && iter->first == pure_id)
        {
            ENSURE(iter->second >= 0);
            return iter->second;
        }
    }
    else if (mixed_id < num_rows())
    {
        // Find unsorted pure ID in the unnormalized row
        auto row_view = row(mixed_id);
        auto iter = std::find_if(row_view.begin(), row_view.end(),
                                 [pure_id](const Matid_Frac &mf)
                                 { return mf.first == pure_id; });
        if (iter != row_view.end())
        {
            ENSURE(iter->second >= 0);
            return iter->second;
        }
    }
    // Not found
    return 0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Concatenate another mix table onto the end of ours
 */
void Mix_Table::extend(const Mix_Table& additional)
{
    REQUIRE(additional.completed() == this->completed());
    REQUIRE(&additional != this);

    // Adjust number of columns (nullop if not completed)
    if (additional.d_num_cols > d_num_cols)
        d_num_cols = additional.d_num_cols;

    // Append table
    d_table.extend(additional.d_table);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a flattened list of all values
 */
Mix_Table::const_View_Frac Mix_Table::values() const
{
    return PROFUGUS_MAKE_STRUCT_VIEW(d_table.values(), second);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a flattened list of all columns
 */
Mix_Table::const_View_Mat Mix_Table::columns() const
{
    return PROFUGUS_MAKE_STRUCT_VIEW(d_table.values(), first);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a flattened list of all row lengths
 */
Mix_Table::const_View_Idx Mix_Table::row_offsets() const
{
    return d_table.offsets();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the COO representation of this table.
 *
 * The input views must be sized to the number of elements in the table.
 */
void Mix_Table::to_coo(View_Mat row, View_Mat col, View_Frac val) const
{
    REQUIRE(row.size() == num_elements());
    REQUIRE(col.size() == num_elements());
    REQUIRE(val.size() == num_elements());

    auto r = row.begin();
    auto c = col.begin();
    auto v = val.begin();

    for (auto row_index : profugus::range(d_table.num_rows()))
    {
        auto view = d_table.row(row_index);
        for (const Matid_Frac& kv : view)
        {
            *r++ = row_index;
            *c++ = kv.first;
            *v++ = kv.second;
        }
    }
    ENSURE(r == row.end());
    ENSURE(c == col.end());
    ENSURE(v == val.end());
}

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
// end of Transcore/material/Mix_Table.cc
//---------------------------------------------------------------------------//
