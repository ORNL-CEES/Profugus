//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Matprop/xs/Mix_Table.i.hh
 * \author Seth R Johnson
 * \date   Thu Feb 04 15:02:21 2016
 * \brief  Mix_Table inline method definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Matprop_xs_Mix_Table_i_hh
#define Matprop_xs_Mix_Table_i_hh

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Un-complete other table on move
 */
Mix_Table::Mix_Table(Mix_Table&& rhs)
    : d_table(std::move(rhs.d_table))
    , d_num_cols(rhs.d_num_cols)
{
    rhs.d_num_cols = 0;
    ENSURE(!rhs.completed());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Un-complete other table on move
 */
Mix_Table& Mix_Table::operator=(Mix_Table&& rhs)
{
    d_table = std::move(rhs.d_table);
    d_num_cols = rhs.d_num_cols;
    rhs.d_num_cols = 0;
    ENSURE(!rhs.completed());
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add a pure matid to this row
 */
void Mix_Table::extend_row(matid_type unmixed_id, value_type fraction)
{
    REQUIRE(fraction >= 0);
    REQUIRE(unmixed_id >= 0);

    if (fraction == 0)
        return;

    d_table.extend_row(std::make_pair(unmixed_id, fraction));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Remove all rows, all columns
 */
void Mix_Table::clear()
{
    d_table.clear();
    d_num_cols = 0;
    ENSURE(!completed());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Swap with another table
 */
void Mix_Table::swap(Mix_Table& rhs)
{
    using std::swap;
    swap(d_table, rhs.d_table);
    swap(d_num_cols, rhs.d_num_cols);
}

//---------------------------------------------------------------------------//
/*!
 * \brief For a given row, access nonzero column/fraction values.
 *
 * If the mix table is not completed, this will return the unnormalized row.
 */
Mix_Table::const_View_Matfrac
Mix_Table::row(matid_type row_id) const
{
    REQUIRE(row_id >= 0 && row_id < num_rows());
    return d_table.row(row_id);
}

//---------------------------------------------------------------------------//
} // end namespace profugus

#endif // Matprop_xs_Mix_Table_i_hh

//---------------------------------------------------------------------------//
// end of Transcore/material/Mix_Table.i.hh
//---------------------------------------------------------------------------//
