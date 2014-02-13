//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/LG_Indexer.hh
 * \author Thomas M. Evans
 * \date   Wednesday February 12 10:30:50 2014
 * \brief  LG_Indexer class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mesh_LG_Indexer_hh
#define mesh_LG_Indexer_hh

#include <vector>

#include "harness/DBC.hh"
#include "utils/Vector_Lite.hh"
#include "utils/Definitions.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class LG_Indexer
 * \brief Convert local indices to global indices on an orthogonal, structured
 * \f$(i,j,k)\f$-indexed mesh or data structure within a set.
 *
 * Using LG_Indexer one can take a locally defined \f$(i,j,k)\f$ and convert
 * it into a global cell index.  Conversely, a globally defined \f$(i,j,k)\f$
 * can be converted into a local cell index/processor pairing, although this
 * is not currently provided.
 *
 * The LG_Indexer is designed to navigate between mesh blocks \e within a
 * set. It is not designed for set-to-set queries. The total number of
 * processor domains is
 * \f[
   N_\mathrm{domains} = N_\mathrm{blocks}\times N_\mathrm{sets}\:.
 * \f]
 */
/*!
 * \example mesh/test/tstLG_Indexer.cc
 *
 * Test of LG_Indexer.
 */
//===========================================================================//

class LG_Indexer
{
  public:
    //@{
    // Useful typedefs.
    typedef def::Vec_Int                      Vec_Int;
    typedef profugus::Vector_Lite<int, 2>     IJ_Set;
    typedef profugus::Vector_Lite<Vec_Int, 2> IJ_Vec;
    typedef std::vector<IJ_Set>               Vec_IJ_Set;
    //@}

  private:
    // >>> DATA

    // Current domain id.
    int d_domain;

    // Processor-block (i,j) indices and block-id.
    IJ_Set d_block;
    int    d_block_id;

    // Number of blocks in (i,j) directions and total blocks.
    IJ_Set d_Nb;
    int    d_num_blocks;

    // Number of cells in each (i,j) block.
    IJ_Vec d_num;

    // Offsets for each (i,j) processor block.
    Vec_Int d_I_offsets;
    Vec_Int d_J_offsets;

    // Global number of cells in (i,j).
    IJ_Set d_num_global;

    // Number of sets and current set-id.
    int d_num_sets;
    int d_set_id;

  public:
    // Constructor.
    LG_Indexer(const Vec_Int &num_I, const Vec_Int &num_J, int Nsets = 1);

    // >>> PUBLIC INTERFACE

    // Set LG_Indexer to a given domain.
    void set_to_domain(int domain);

    // >>> PUBLIC INTERFACE
    // Return the current domain.
    int current_domain() const { return d_domain; }

    //! Number of cells in \f$(i,j)\f$ on current domain.
    int num_cells(int d) const { return d_num[d][d_block[d]]; }

    // Convert local (i,j,k) to global cell index.
    inline int l2g(int i, int j, int k) const;

    // Convert local (i,j,k) to local cell index.
    inline int l2l(int i, int j, int k) const;

    // Convert local cell index to local (i,j,k).
    inline void l2l(int cell, int &i, int &j, int &k) const;

    // Convert global (i,j,k) to a global cell index.
    inline int g2g(int i, int j, int k) const;

    // Convert global cell to global (i,j,k).
    inline void g2g(int cell, int &i, int &j, int &k) const;

    // Convert (i,j) from local to global cell indices.
    inline IJ_Set convert_to_global(int i, int j) const;

    // Convert (i,j) from global to local indices.
    inline IJ_Set convert_to_local(int i, int j) const;

    // Return the set index.
    int set() const { return d_set_id; }

    // Return the number of sets.
    int num_sets() const { return d_num_sets; }

    // Return the block index.
    int block() const { return d_block_id; }

    //! Number of blocks (domains) in a set.
    int num_blocks() const { return d_num_blocks; }

    //! Number of blocks (domains) in \f$(i,j)\f$ directions.
    int num_blocks(int d) const { Require (d < 3); return d_Nb[d]; }

    //! Number of global cells in \f$(i,j)\f$.
    int num_global(int d) const { return d_num_global[d]; }

    //! Convert a block and set to a domain
    inline int domain(int block, int set) const;

    // Number of cells in \f$(i,j)\f$ directions.
    const Vec_Int& num_cells_per_block(int d) const { return d_num[d]; }

    // Return the \f$(i,j)\f$ offsets for the block on the current node.
    inline int offset(int dir) const;

    // Obtain the domains and on-processor cells that span a global range
    void global_to_manylocal(const IJ_Set& gbegin, const IJ_Set& gend,
                             int set, Vec_Int& domains, Vec_IJ_Set& lbegins,
                             Vec_IJ_Set& lends) const;

  private:
    // >>> DATA AND IMPLEMENTATION

    // Total number of domains (processors).
    int d_nodes;
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Convert a local \f$(i,j,k)\f$ to a global cell index.
 *
 * \return global cell index in the range [0, N)
 */
int LG_Indexer::l2g(int i,
                    int j,
                    int k) const
{
    using def::I;
    using def::J;

    Require (i >= 0);
    Require (j >= 0);
    Require (k >= 0);
    Require (i < d_num[I][d_block[I]]);
    Require (j < d_num[J][d_block[J]]);

    return (i + d_I_offsets[d_block[I]]) + d_num_global[I] * (
        j + d_J_offsets[d_block[J]] + k * d_num_global[J]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert a local \f$(i,j,k)\f$ to a local cell index.
 *
 * \return local cell index in the range [0, N)
 */
int LG_Indexer::l2l(int i,
                    int j,
                    int k) const
{
    using def::I;
    using def::J;

    Require (i >= 0);
    Require (j >= 0);
    Require (k >= 0);
    Require (i < d_num[I][d_block[I]]);
    Require (j < d_num[J][d_block[J]]);

    return (i + j * d_num[I][d_block[I]] + k * d_num[I][d_block[I]] *
            d_num[J][d_block[J]]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert a local cell index to a local \f$(i,j,k)\f$.
 */
void LG_Indexer::l2l(int  cell,
                     int &i,
                     int &j,
                     int &k) const
{
    using def::I;
    using def::J;

    // store NxNy
    int NxNy = d_num[I][d_block[I]] * d_num[J][d_block[J]];

    // get k
    k = cell / NxNy;

    // get j
    j = (cell - k * NxNy) / d_num[I][d_block[I]];

    // get i
    i = cell - j * d_num[I][d_block[I]] - k * NxNy;

    Ensure (cell == l2l(i, j, k));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert a global \f$(i,j,k)\f$ to a global cell index.
 *
 * \return global cell index in the range [0, N)
 */
int LG_Indexer::g2g(int i,
                    int j,
                    int k) const
{
    using def::I;
    using def::J;

    Require (i >= 0);
    Require (j >= 0);
    Require (k >= 0);
    Require (i < d_num_global[I]);
    Require (j < d_num_global[J]);

    return i + d_num_global[I] * (j + k * d_num_global[J]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert a global cell index to a global \f$(i,j,k)\f$.
 */
void LG_Indexer::g2g(int  cell,
                     int &i,
                     int &j,
                     int &k) const
{
    using def::I;
    using def::J;

    // store NxNy
    int NxNy = d_num_global[I] * d_num_global[J];

    // get k
    k = cell / NxNy;

    // get j
    j = (cell - k * NxNy) / d_num_global[I];

    // get i
    i = cell - j * d_num_global[I] - k * NxNy;

    Ensure (cell == g2g(i, j, k));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert \f$(i,j)\f$ from local to global.
 *
 * \return IJ_Set of global \f$(i,j)\f$ indices
 */
LG_Indexer::IJ_Set LG_Indexer::convert_to_global(int i,
                                                 int j) const
{
    using def::I;
    using def::J;

    Require (i >= 0);
    Require (j >= 0);
    Require (i < d_num[I][d_block[I]]);
    Require (j < d_num[J][d_block[J]]);

    return IJ_Set(i + d_I_offsets[d_block[I]], j + d_J_offsets[d_block[J]]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert \f$(i,j)\f$ from global to local.
 *
 * \return IJ_Set of global \f$(i,j)\f$ indices, entries will be -1 if the
 * global indices are not local to the set node
 */
LG_Indexer::IJ_Set LG_Indexer::convert_to_local(int i,
                                                int j) const
{
    using def::I;
    using def::J;

    Require (i >= 0);
    Require (j >= 0);
    Require (i < d_num_global[I]);
    Require (j < d_num_global[J]);

    // default return value
    IJ_Set ij(-1, -1);

    // both indices must be local to return a value, check I first
    if (i >= d_I_offsets[d_block[I]] &&
        i < d_I_offsets[d_block[I]] + d_num[I][d_block[I]])
    {
        if (j >= d_J_offsets[d_block[J]] &&
            j < d_J_offsets[d_block[J]] + d_num[J][d_block[J]])
        {
            ij[I] = i - d_I_offsets[d_block[I]];
            ij[J] = j - d_J_offsets[d_block[J]];

            Ensure (ij[I] >= 0 && ij[I] < num_cells(I));
            Ensure (ij[J] >= 0 && ij[J] < num_cells(J));
        }
    }

    return ij;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the domain that a given block and set belong to.
 */
int LG_Indexer::domain(int block, int set) const
{
    Require (block < num_blocks());
    Require (set < num_sets());

    return set*num_blocks() + block;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the \f$(i,j)\f$ offsets for the block on the current domain.
 *
 * The offsets are the number of cells in all blocks leading up to the current
 * block such that:
 * \f[
   \mbox{global\_cell} = (i + \mbox{offset}_i) + (j + \mbox{offset}_j)\times
   N_i + k\times N_i\times N_j\:,
 * \f]
 * where \f$(i,j,k)\f$ are defined locally on a block over the range \f$(0,
 * N]\f$, and \f$(N_i, N_j)\f$ are the number of cells in the block in the \f$
 * i\f$ and \f$ j\f$ directions, respectively.
 */
int LG_Indexer::offset(int dir) const
{
    Require (dir == def::I || dir == def::J);
    if (dir == def::I) return d_I_offsets[d_block[def::I]];
    return d_J_offsets[d_block[def::J]];
}

} // end namespace profugus

#endif // mesh_LG_Indexer_hh

//---------------------------------------------------------------------------//
//              end of mesh/LG_Indexer.hh
//---------------------------------------------------------------------------//
