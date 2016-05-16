//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/utils/Regular_Indexer.hh
 * \author Seth R Johnson
 * \date   Mon Aug 10 16:30:18 2015
 * \brief  Regular_Indexer class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Regular_Indexer_hh
#define Utils_utils_Regular_Indexer_hh

#include <cstddef>
#include "Vector_Lite.hh"

namespace profugus
{
//===========================================================================//
/*!
 * \class Regular_Indexer
 * \brief Indexing implementation for hyperslab views.
 *
 * This handles the per-element strides associated with the dimensioning.
 * (Additional striding, e.g. over reactions, mean/variance, unknowns, can be
 * provided implicitly by the stride in the view field iterator.)
 *
 * The dimensions are stored in row-major (C-like) order, with the
 * fastest-moving dimension on the right (index.back())
 *
 * \warning Because the dimensions are row-major, the index (like a Dim_Vector)
 * will be ordered in opposite order of a Mesh_KBA Dim_Vector, which uses the
 * Fortran-like column-major (I,J,K)!! This is also the case for the
 * profugus::Regular_Grid class, which follows the KBA indexing.
 */
/*!
 * \example utils/test/tstHyperslab_Indexer.cc
 *
 * Test of tstHyperslab_Indexer.
 */
//===========================================================================//

template<typename T, std::size_t D>
class Regular_Indexer
{
  public:
    //@{
    //! Typedefs
    using size_type     = T;
    using index_type    = Vector_Lite<size_type, D>;
    using subindex_type = Vector_Lite<size_type, D-1>;
    //@}

  private:
    // >>> DATA

    //! Size of each dimension
    index_type d_dims;

  public:
    // >>> CONSTRUCTOR

    // Construct with empty dimensions
    inline Regular_Indexer();

    // Construct with dimensionality
    explicit inline Regular_Indexer(index_type dims);

    // >>> ACCESSORS

    //! Access dimensions
    const index_type& dims() const { return d_dims; }

    // Access the total number of elements
    inline size_type size() const;

    // Get dimensions for a slice along the major index
    inline subindex_type major_slice_dims() const;

    // Get dimensions for a slice along the minor index
    inline subindex_type minor_slice_dims() const;

    // Get the index in the flattened array of the element at i
    inline size_type index(const index_type &i) const;

    // Get the index in the flattened array of the element at i
    inline index_type index(size_type i) const;

    // >>> SERIALIZE

    //! Serialization
    template<class Archiver>
    void serialize(Archiver &ar)
    {
        ar & d_dims;
    }
};

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "Regular_Indexer.i.hh"
//---------------------------------------------------------------------------//
#endif // Utils_utils_Regular_Indexer_hh

//---------------------------------------------------------------------------//
// end of Regular_Indexer.hh
//---------------------------------------------------------------------------//
