//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/utils/Regular_Indexer_Impl.hh
 * \author Seth R Johnson
 * \date   Mon Nov 23 10:47:20 2015
 * \brief  Regular_Indexer_Impl class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Regular_Indexer_Impl_hh
#define Utils_utils_Regular_Indexer_Impl_hh

#include <cstddef>
#include "Utils/harness/DBC.hh"
#include "Vector_Lite.hh"

namespace profugus
{
namespace detail
{

//===========================================================================//
/*!
 * \struct Regular_Indexer_Impl
 *
 * Implementation of indexing to allow for specicalizations.
 *
 * In GCC5 with -O2, these automatically unroll up to and including 3-dimension
 * vectors.
 */
//===========================================================================//

template<typename T, std::size_t D>
struct Regular_Indexer_Impl
{
    using size_type  = T;
    using index_type = ::profugus::Vector_Lite<size_type, D>;

    //! Convert index vector to flattened index
    static size_type index(const index_type& dims, const index_type& i)
    {
        size_type result = i[0];

        for (std::size_t d = 1u; d < D; ++d)
        {
            result = dims[d] * result + i[d];
        }

        return result;
    }

    /*! Convert flattened index to index vector.
     *
     * This method is O(D^2), but it has a lower constant since it avoids
     * additional integer divisions.
     */
    static index_type index(const index_type& dims, size_type i)
    {
        index_type result;

        for (std::size_t d = 0u; d < D - 1u; ++d)
        {
            size_type stride = 1u;
            for (std::size_t dp = d + 1u; dp < D; ++dp)
            {
                stride *= dims[dp];
            }

            const size_type temp = i / stride;
            result[d] = temp;
            i -= temp * stride;
        }

        result[D - 1u] = i;
        return result;
    }
};

//===========================================================================//
// 0-D specialization: no index methods!
//===========================================================================//

template<typename T>
struct Regular_Indexer_Impl<T, 0u>
{
    using size_type = T;

    template<class I>
    static size_type index(const I&, const I&)
    {
        static_assert(sizeof(T) == 0, "Can't instantiate 0D indexer!");
    }

    template<class I>
    static I index(const I&, size_type)
    {
        static_assert(sizeof(T) == 0, "Can't instantiate 0D indexer!");
    }
};

//---------------------------------------------------------------------------//
} // end namespace profugus::detail
} // end namespace profugus

//---------------------------------------------------------------------------//
#endif // Utils_utils_Regular_Indexer_Impl_hh

//---------------------------------------------------------------------------//
// end of Regular_Indexer_Impl.hh
//---------------------------------------------------------------------------//
