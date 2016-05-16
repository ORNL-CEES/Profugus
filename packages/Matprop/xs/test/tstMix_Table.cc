//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Matprop/xs/test/tstMix_Table.cc
 * \author Seth R Johnson
 * \date   Thu Feb 04 18:38:16 2016
 * \brief  Tests for class Mix_Table.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Mix_Table.hh"

#include "Utils/gtest/utils_gtest.hh"
#include "Utils/utils/View_Field.hh"

using profugus::Mix_Table;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class MixTableTest : public ::testing::Test
{
  protected:
    using matid_type = Mix_Table::matid_type;
    using value_type = Mix_Table::value_type;

    void expand(const Mix_Table& m, matid_type row)
    {
        ASSERT_TRUE(m.completed());
        ASSERT_LT(row, m.num_rows());
        ids.clear();
        fracs.clear();
        for (const auto& id_frac : m.row(row))
        {
            ids.push_back(id_frac.first);
            fracs.push_back(id_frac.second);
        }
    }

  protected:
    // >>> DATA
    std::vector<matid_type> ids;
    std::vector<value_type> fracs;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(MixTableTest, construction)
{
    Mix_Table table;
    EXPECT_FALSE(table.completed());

    table.start_row();
    table.extend_row(1, 1.0);
    table.extend_row(3, 3.0);
    table.start_row();
    table.extend_row(2, 1.0);
    table.extend_row(3, 1.0);
    table.start_row();
    // Out-of-order insertion
    table.extend_row(5, 0.5);
    table.extend_row(3, 1.0);
    table.extend_row(0, 0.5);
    table.start_row();

    EXPECT_SOFT_EQ(0.0, table.find_fraction(10, 0));
    EXPECT_SOFT_EQ(0.0, table.find_fraction(0, 10));
    EXPECT_SOFT_EQ(1.0, table.find_fraction(0, 1));
    EXPECT_SOFT_EQ(0.5, table.find_fraction(2, 0));

    table.complete();
    EXPECT_EQ(4, table.num_rows());
    EXPECT_EQ(6, table.num_cols());

    expand(table, 0);
    {
        const matid_type expected_ids[] = {1, 3};
        const value_type expected_fracs[] = {0.25, 0.75};
        EXPECT_VEC_EQ(expected_ids, ids);
        EXPECT_VEC_SOFT_EQ(expected_fracs, fracs);
    }

    expand(table, 1);
    {
        const matid_type expected_ids[] = {2, 3};
        const value_type expected_fracs[] = {0.5, 0.5};
        EXPECT_VEC_EQ(expected_ids, ids);
        EXPECT_VEC_SOFT_EQ(expected_fracs, fracs);
    }

    expand(table, 2);
    {
        const matid_type expected_ids[] = {0, 3, 5};
        const value_type expected_fracs[] = {0.25, 0.5, 0.25};
        EXPECT_VEC_EQ(expected_ids, ids);
        EXPECT_VEC_SOFT_EQ(expected_fracs, fracs);
    }

    expand(table, 3);
    {
        EXPECT_EQ(0, ids.size());
        EXPECT_EQ(0, fracs.size());
    }

#ifdef REQUIRE_ON
    EXPECT_THROW(table.find_fraction(10, 0), profugus::assertion);
    EXPECT_THROW(table.find_fraction(0, 10), profugus::assertion);
#endif
    EXPECT_SOFT_EQ(0.25, table.find_fraction(0, 1));
    EXPECT_SOFT_EQ(0.25, table.find_fraction(2, 0));
    EXPECT_SOFT_EQ(0.0, table.find_fraction(1, 0));
    EXPECT_SOFT_EQ(0.5, table.find_fraction(2, 3));
}

TEST_F(MixTableTest, from_coo)
{
    // Note: last row is unnormalized
    int rows[]    = {0,   1,  1,  1, 2, 4};
    int cols[]    = {1,   0,  1,  2, 2, 3};
    double vals[] = {1., .1, .7, .2, 1., 2.};

    Mix_Table table = Mix_Table::from_coo(profugus::make_view(rows),
                                          profugus::make_view(cols),
                                          profugus::make_view(vals));
    EXPECT_TRUE(table.completed());

    const double expected_values[] = {1, 0.1, 0.7, 0.2, 1, 1, };
    const int expected_columns[] = {1, 0, 1, 2, 2, 3, };
    const std::size_t expected_row_offsets[] = {0, 1, 4, 5, 5, 6, };

    EXPECT_VEC_SOFT_EQ(expected_values, table.values());
    EXPECT_VEC_EQ(expected_columns, table.columns());
    EXPECT_VEC_EQ(expected_row_offsets, table.row_offsets());

    auto num_els = table.num_elements();
    std::vector<int> actual_rows(num_els, -1);
    std::vector<int> actual_cols(num_els, -1);
    std::vector<double> actual_vals(num_els, -1);

    // Renormalize expected value
    vals[5] = 1.0;
    table.to_coo(profugus::make_view(actual_rows),
                 profugus::make_view(actual_cols),
                 profugus::make_view(actual_vals));

    EXPECT_VEC_EQ(rows, actual_rows);
    EXPECT_VEC_EQ(cols, actual_cols);
    EXPECT_VEC_SOFT_EQ(vals, actual_vals);

    //PRINT_EXPECTED(table.values());
    //PRINT_EXPECTED(table.columns());
    //PRINT_EXPECTED(table.row_offsets());
}

//---------------------------------------------------------------------------//
// end of Transcore/material/test/tstMix_Table.cc
//---------------------------------------------------------------------------//
