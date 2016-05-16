//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/utils/test/tstRegular_Indexer.cc
 * \author Seth R Johnson
 * \date   Mon Aug 10 16:33:10 2015
 * \brief  Regular_Indexer class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Regular_Indexer.hh"

#include "Utils/gtest/utils_gtest.hh"

using profugus::Regular_Indexer;
using std::size_t;

//---------------------------------------------------------------------------//
// HYPERSLAB_INDEXER TESTS
//---------------------------------------------------------------------------//
TEST(HyperslabIndexerTest, oned)
{
    typedef Regular_Indexer<size_t, 1>  Indexer_t;
    typedef Indexer_t::index_type index_type;

    index_type dims(3);

    Indexer_t indexer(dims);

    EXPECT_VEC_EQ(dims, indexer.dims());
    EXPECT_EQ(3, indexer.size());

    index_type i;
    i = index_type(0); EXPECT_EQ(0, indexer.index(i));
    i = index_type(1); EXPECT_EQ(1, indexer.index(i));
    i = index_type(2); EXPECT_EQ(2, indexer.index(i));
#ifdef REQUIRE_ON
    i = index_type(4); EXPECT_THROW(indexer.index(i), profugus::assertion);
#endif
}

//---------------------------------------------------------------------------//
TEST(HyperslabIndexerTest, twod)
{
    typedef Regular_Indexer<int, 2> Indexer_t;
    typedef Indexer_t::index_type     index_type;
    typedef Indexer_t::subindex_type  subindex_type;

    index_type dims(2, 3);

    Indexer_t indexer(dims);

    EXPECT_VEC_EQ(dims, indexer.dims());
    EXPECT_EQ(2 * 3, indexer.size());
    EXPECT_EQ(subindex_type(3), indexer.major_slice_dims());

    int ctr = 0;
    for (int j = 0; j != 2; ++j)
    {
        for (int i = 0; i != 3; ++i, ++ctr)
        {
            EXPECT_EQ(ctr, indexer.index(index_type(j,i)))
                << "Failure with (j,i) = (" << j << "," << i << ")";
        }
    }
#ifdef REQUIRE_ON
    EXPECT_THROW(indexer.index(index_type(3,2)), profugus::assertion);
    EXPECT_THROW(indexer.index(index_type(2,4)), profugus::assertion);
#endif
}

//---------------------------------------------------------------------------//
TEST(HyperslabIndexerTest, threed)
{
    typedef Regular_Indexer<unsigned int, 3> Indexer_t;
    typedef Indexer_t::index_type              index_type;
    typedef Indexer_t::subindex_type           subindex_type;

    index_type dims(2, 3, 4);

    Indexer_t indexer(dims);

    EXPECT_VEC_EQ(dims, indexer.dims());
    EXPECT_EQ(2 * 3 * 4, indexer.size());
    EXPECT_EQ(subindex_type(3, 4), indexer.major_slice_dims());

    int ctr = 0;
    for (int k = 0; k != 2; ++k)
    {
        for (int j = 0; j != 3; ++j)
        {
            for (int i = 0; i != 4; ++i, ++ctr)
            {
                EXPECT_EQ(ctr, indexer.index(index_type(k,j,i)))
                    << "Failure with (k,j,i) = ("
                    << k << "," << j << "," << i << ")";
            }
        }
    }
#ifdef ENSURE_ON
    EXPECT_THROW(indexer.index(index_type(3,0,1)), profugus::assertion);
#endif
#ifdef REQUIRE_ON
    EXPECT_THROW(indexer.index(index_type(0,5,2)), profugus::assertion);
    EXPECT_THROW(indexer.index(index_type(3,2,4)), profugus::assertion);
#endif
}

//---------------------------------------------------------------------------//
// end of tstRegular_Indexer.cc
//---------------------------------------------------------------------------//
