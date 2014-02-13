//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/LG_Indexer.cc
 * \author Thomas M. Evans
 * \date   Tue Aug 28 14:43:58 2007
 * \brief  LG_Indexer member definitions.
 * \note   Copyright (C) 2007 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include "comm/global.hh"
#include "LG_Indexer.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * Here, \c Nb(I) is the number of processor-blocks in the \e I -direction, \c
 * Nb(J) is the number of processors in the \e J -direction.
 *
 * The LG_Indexer is constructed to work on the current domain.  It can be set
 * to work on any node using set_to_domain().
 *
 * \param num_I vector of number of cells in \e I -direction per block;
 * defined over the range \c [0,Nb(I))
 *
 * \param num_J vector of number of cells in \e J -direction per block;
 * defined over the range \c [0,Nb(J))
 *
 * \param Nsets number of sets (optional, defaults to 1)
 */
LG_Indexer::LG_Indexer(const Vec_Int &num_I,
                       const Vec_Int &num_J,
                       int            Nsets)
    : d_Nb(num_I.size(), num_J.size())
    , d_num_blocks(d_Nb[def::I] * d_Nb[def::J])
    , d_num(num_I, num_J)
    , d_I_offsets(num_I.size() + 1, 0)
    , d_J_offsets(num_J.size() + 1, 0)
    , d_num_sets(Nsets)
    , d_nodes(profugus::nodes())
{
    using def::I;
    using def::J;

    Require (d_num[I].size() * d_num[J].size() == d_num_blocks);
    Require (d_num_blocks * d_num_sets == d_nodes);

    // set LG_Indexer to current domain (processor node)
    set_to_domain(profugus::node());

    // calculate the offsets
    for (int i = 1; i <= d_Nb[I]; i++)
    {
        d_I_offsets[i] = d_I_offsets[i-1] + d_num[I][i-1];
    }
    for (int j = 1; j <= d_Nb[J]; j++)
    {
        d_J_offsets[j] = d_J_offsets[j-1] + d_num[J][j-1];
    }

    // set the global number of cells in (i,j)
    d_num_global[I] = d_I_offsets.back();
    d_num_global[J] = d_J_offsets.back();

    Ensure (d_num_global[I] ==
            std::accumulate(d_num[I].begin(), d_num[I].end(), 0));
    Ensure (d_num_global[J] ==
            std::accumulate(d_num[J].begin(), d_num[J].end(), 0));
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set the LG_Indexer to provide indices for a given domain.
 */
void LG_Indexer::set_to_domain(int domain)
{
    Require (domain >= 0 && domain < d_nodes);

    using std::accumulate;
    using def::I;
    using def::J;
    using def::K;

    // assign domain
    d_domain = domain;

    // calculate the set id
    d_set_id = d_domain / (d_num_blocks);
    Check (d_set_id >= 0 && d_set_id < d_num_sets);

    // calculate the block id
    d_block_id = d_domain - d_set_id * d_num_blocks;
    Check (d_block_id >= 0 && d_block_id < d_num_blocks);

    // calculate the processor index
    d_block[J] = d_block_id / d_Nb[I];
    d_block[I] = d_block_id - d_block[J] * d_Nb[I];

    Ensure (d_block[J] < d_Nb[J]);
    Ensure (d_block[I] < d_Nb[I]);
    Ensure (d_block[I] + d_block[J] * d_Nb[I] == d_block_id);
    Ensure (d_block_id + d_set_id * d_num_blocks == d_domain);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Obtain the domains and on-processor cells that span a global range
 *
 * This function will determine which domains span the given global cells for
 * the given set.
 *
 * The input range must contain more than one cell.
 */
void LG_Indexer::global_to_manylocal(const IJ_Set& gbegin,
                                     const IJ_Set& gend,
                                     int           set,
                                     Vec_Int&      domains,
                                     Vec_IJ_Set&   lbegins,
                                     Vec_IJ_Set&   lends) const
{
    using def::I; using def::J;

    Require(0         <= gbegin[I]);     Require(0         <= gbegin[J]);
    Require(gbegin[I] <  gend[I]);       Require(gbegin[J] <  gend[J]);
    Require(gend[I]   <= num_global(I)); Require(gend[J]   <= num_global(J));
    Require(set < num_sets());

    domains.clear();
    lbegins.clear();
    lends.clear();

    // find begin and end block indices in i and j
    int ib_begin = std::lower_bound(
            d_I_offsets.begin(),
            d_I_offsets.end(),
            gbegin[def::I]) - d_I_offsets.begin();
    if (d_I_offsets[ib_begin] != gbegin[def::I])
        --ib_begin;
    Check(0 <= ib_begin && ib_begin < num_blocks(I));

    int ib_end = std::lower_bound(
            d_I_offsets.begin(),
            d_I_offsets.end(),
            gend[def::I]) - d_I_offsets.begin();
    Check(0 < ib_end && ib_end <= num_blocks(I));
    Check(ib_begin < ib_end);

    int jb_begin = std::lower_bound(
            d_J_offsets.begin(),
            d_J_offsets.end(),
            gbegin[def::J]) - d_J_offsets.begin();
    if (d_J_offsets[jb_begin] != gbegin[def::J])
        --jb_begin;
    Check(0 <= jb_begin && jb_begin < num_blocks(J));

    int jb_end = std::lower_bound(
            d_J_offsets.begin(),
            d_J_offsets.end(),
            gend[def::J]) - d_J_offsets.begin();
    Check(0 < jb_end && jb_end <= num_blocks(J));
    Check(jb_begin < jb_end);

    const unsigned int num_domains = (jb_end - jb_begin) * (ib_end - ib_begin);
    domains.reserve(num_domains);
    lbegins.reserve(num_domains);
    lends.reserve(num_domains);

    IJ_Set lbegin;
    IJ_Set lend;

    for (int jb = jb_begin; jb < jb_end; ++jb)
    {
        // Set local beginning and ending indices for J
        lbegin[J] = 0;
        if (d_J_offsets[jb] < gbegin[J])
            lbegin[J] = gbegin[J] - d_J_offsets[jb];

        lend[J] = d_num[J][jb];
        if (gend[J] < d_J_offsets[jb + 1])
            lend[J] = gend[J] - d_J_offsets[jb];

        Check(0 <= lbegin[J]);
        Check(lbegin[J] < lend[J]);
        Check(lend[J] <= d_num[J][jb]);

        for (int ib = ib_begin; ib < ib_end; ++ib)
        {
            // Set local beginning and ending indices for I
            lbegin[I] = 0;
            if (d_I_offsets[ib] < gbegin[I])
                lbegin[I] = gbegin[I] - d_I_offsets[ib];

            lend[I] = d_num[I][ib];
            if (gend[I] < d_I_offsets[ib + 1])
                lend[I] = gend[I] - d_I_offsets[ib];

            Check(0 <= lbegin[I]);
            Check(lbegin[I] < lend[I]);
            Check(lend[I] <= d_num[I][ib]);

            // Add the domain
            domains.push_back(domain(ib + jb * num_blocks(I), set));
            // Add the start/end indices
            lbegins.push_back(lbegin);
            lends.push_back(lend);
        }
    }

    Ensure(domains.size() == num_domains);
    Ensure(domains.size() == lbegins.size());
    Ensure(domains.size() == lends.size());
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//              end of LG_Indexer.cc
//---------------------------------------------------------------------------//
