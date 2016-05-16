//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/utils/Regular_Indexer.i.hh
 * \author Seth R Johnson
 * \date   Mon Aug 10 16:30:18 2015
 * \brief  Regular_Indexer inline method definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Regular_Indexer_i_hh
#define Utils_utils_Regular_Indexer_i_hh

#include <type_traits>
#include "Regular_Indexer_Impl.hh"

namespace profugus
{
//---------------------------------------------------------------------------//
// HYPERSLAB_INDEXER
//---------------------------------------------------------------------------//
/*!
 * \brief Construct with all dimensions set to zero.
 *
 * This is used for the default constructor of a hyperslab view.
 */
template<typename T, std::size_t D>
Regular_Indexer<T,D>::Regular_Indexer()
    : d_dims(0)
{
    /* * */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct with dimensions.
 */
template<typename T, std::size_t D>
Regular_Indexer<T,D>::Regular_Indexer(index_type dims)
    : d_dims(dims)
{
    static_assert(D > 0, "Regular_Indexer dimension must be positive");
    static_assert(std::is_integral<T>::value,
                  "Regular_Indexer must be templated on an integral type");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the total number of elements.
 */
template<typename T, std::size_t D>
T Regular_Indexer<T,D>::size() const
{
    return std::accumulate(d_dims.begin(), d_dims.end(),
                           static_cast<size_type>(1),
                           std::multiplies<size_type>());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the dimensions if the major axis is sliced off.
 */
template<typename T, std::size_t D>
typename Regular_Indexer<T,D>::subindex_type
Regular_Indexer<T,D>::major_slice_dims() const
{
    subindex_type subindex;
    for (std::size_t d = 1u; d != D; ++d)
    {
        subindex[d-1] = d_dims[d];
    }

    return subindex;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the dimensions if the minor axis is sliced off.
 */
template<typename T, std::size_t D>
typename Regular_Indexer<T,D>::subindex_type
Regular_Indexer<T,D>::minor_slice_dims() const
{
    subindex_type subindex;
    for (std::size_t d = 0u; d != D - 1u; ++d)
    {
        subindex[d] = d_dims[d];
    }

    return subindex;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the index in the flattened array of the element at i
 */
template<typename T, std::size_t D>
T Regular_Indexer<T,D>::index(const index_type &i) const
{
#ifdef REQUIRE_ON
    for (std::size_t d = 0u; d != D; ++d)
        REQUIRE(i[d] <= d_dims[d]);
#endif

    using Indexer_Impl_t = ::profugus::detail::Regular_Indexer_Impl<T, D>;
    auto result = Indexer_Impl_t::index(d_dims, i);
    ENSURE(result <= size());
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the cardinal index of flattened element i
 */
template<typename T, std::size_t D>
typename Regular_Indexer<T,D>::index_type
Regular_Indexer<T,D>::index(size_type i) const
{
    using Indexer_Impl_t = ::profugus::detail::Regular_Indexer_Impl<T, D>;
    return Indexer_Impl_t::index(d_dims, i);
}

//---------------------------------------------------------------------------//
} // end namespace profugus

#endif // Utils_utils_Regular_Indexer_i_hh

//---------------------------------------------------------------------------//
// end of Regular_Indexer.i.hh
//---------------------------------------------------------------------------//
