//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/Hyperslab_View.i.hh
 * \author Seth R Johnson
 * \date   Mon Dec 08 15:25:35 2014
 * \brief  Member definitions of class Hyperslab_View.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Hyperslab_View_i_hh
#define Utils_utils_Hyperslab_View_i_hh

#include <functional>

namespace profugus
{
//---------------------------------------------------------------------------//
// HYPERSLAB_VIEW
//---------------------------------------------------------------------------//
/*!
 * \brief Construct with empty view
 */
template<class T, size_t D>
Hyperslab_View<T,D>::Hyperslab_View()
    : Base()
    , d_indexer()
{
    ENSURE(size() == 0);
    ENSURE(size() == d_indexer.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct with data pointers, dimensionality.
 */
template<class T, size_t D>
Hyperslab_View<T,D>::Hyperslab_View(
        view data,
        index_type dims)
    : Base(data)
    , d_indexer(dims)
{
    REQUIRE(data.size() == d_indexer.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a single element of the container.
 */
template<class T, size_t D>
typename Hyperslab_View<T,D>::reference
Hyperslab_View<T,D>::operator[] (const index_type& index)
{
    size_type i = d_indexer.index(index);
    REQUIRE(i < size());
    return *(begin() + i);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a single element of the container.
 */
template<class T, size_t D>
typename Hyperslab_View<T,D>::const_reference
Hyperslab_View<T,D>::operator[] (const index_type& index) const
{
    size_type i = d_indexer.index(index);
    REQUIRE(i < size());
    return *(begin() + i);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a slice of the hyperslab (most major index).
 *
 * \warning There will be several integer multiplies and some on-stack integer
 * copies associated with this operation. Save and reuse these slices whenever
 * possible.
 */
template<class T, size_t D>
typename Hyperslab_View<T,D>::subview_type
Hyperslab_View<T,D>::major_slice(size_type index)
{
    REQUIRE(index < d_indexer.dims().front());

    // Create slice dimensions
    const auto &subindex = d_indexer.major_slice_dims();

    // Slice the data
    size_type width = std::accumulate(subindex.begin(), subindex.end(), 1u,
                                      std::multiplies<size_type>());
    view subview(data() + width * index,
                 data() + width * (index + 1));

    return subview_type(subview, subindex);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a slice of the hyperslab (most major index).
 */
template<class T, size_t D>
typename Hyperslab_View<T,D>::const_subview_type
Hyperslab_View<T,D>::major_slice(size_type index) const
{
    REQUIRE(index < d_indexer.dims().front());

    // Create slice dimensions
    const auto &subindex = d_indexer.major_slice_dims();

    // Slice the data
    size_type width = std::accumulate(subindex.begin(), subindex.end(), 1u,
                                      std::multiplies<size_type>());
    const_view subview(begin().ptr() + width * index,
                       begin().ptr() + width * (index + 1));

    return const_subview_type(subview, subindex);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a slice of the hyperslab (most minor index).
 *
 * \warning This returns a strided view field; it will not be contiguous in
 * memory.
 *
 * \warning There will be several integer multiplies and some on-stack integer
 * copies associated with this operation. Save and reuse these slices whenever
 * possible.
 */
template<class T, size_t D>
typename Hyperslab_View<T,D>::subview_type
Hyperslab_View<T,D>::minor_slice(size_type index)
{
    REQUIRE(index < d_indexer.dims().back());

    // Create slice dimensions
    const auto &subindex = d_indexer.minor_slice_dims();

    // Slice the data
    view subview = strided_slice(index,
                                 index + size(),
                                 d_indexer.dims().back());

    return subview_type(subview, subindex);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a slice of the hyperslab (most minor index).
 *
 * \warning This returns a strided view field; it will not be contiguous in
 * memory.
 *
 * \warning There will be several integer multiplies and some on-stack integer
 * copies associated with this operation. Save and reuse these slices whenever
 * possible.
 */
template<class T, size_t D>
typename Hyperslab_View<T,D>::const_subview_type
Hyperslab_View<T,D>::minor_slice(size_type index) const
{
    REQUIRE(index < d_indexer.dims().back());

    // Create slice dimensions
    const auto &subindex = d_indexer.minor_slice_dims();

    // Slice the data
    const_view subview = strided_slice(index,
                                       index + size(),
                                       d_indexer.dims().back());

    return const_subview_type(subview, subindex);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator
 */
template<class T, size_t D>
Hyperslab_View<T,D>&
Hyperslab_View<T,D>::operator=(const Hyperslab_View<T,D>& rhs)
{
    Base::operator=(rhs);
    d_indexer = rhs.d_indexer;
    ENSURE(size() == d_indexer.size());
    return *this;
}

//---------------------------------------------------------------------------//
// CONST_HYPERSLAB_VIEW
//---------------------------------------------------------------------------//
/*!
 * \brief Construct with empty view
 */
template<class T, size_t D>
const_Hyperslab_View<T,D>::const_Hyperslab_View()
    : Base()
    , d_indexer()
{
    ENSURE(size() == 0);
    ENSURE(size() == d_indexer.size());
}
//---------------------------------------------------------------------------//
/*!
 * \brief Construct with data pointers, dimensionality.
 *
 * The stride of the data gets put into the "minor" index
 */
template<class T, size_t D>
const_Hyperslab_View<T,D>::const_Hyperslab_View(
        const_view data,
        index_type dims)
    : Base(data)
    , d_indexer(dims)
{
    REQUIRE(size() == d_indexer.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct implicitly from non-const view.
 *
 * The stride of the data gets put into the "minor" index
 */
template<class T, size_t D>
const_Hyperslab_View<T,D>::const_Hyperslab_View(
        const Hyperslab_View<T,D> &rhs)
    : Base(static_cast<const View_Field<T>&>(rhs))
    , d_indexer(rhs.d_indexer)
{
    ENSURE(stride() == rhs.stride());
    ENSURE(size() == rhs.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a single element of the container.
 */
template<class T, size_t D>
typename const_Hyperslab_View<T,D>::const_reference
const_Hyperslab_View<T,D>::operator[] (const index_type& index) const
{
    size_type i = d_indexer.index(index);
    REQUIRE(i < size());
    return *(begin() + i);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a slice of the hyperslab (most major index).
 *
 * \warning There will be several integer multiplies and some on-stack integer
 * copies associated with this operation. Save and reuse these slices whenever
 * possible.
 */
template<class T, size_t D>
typename const_Hyperslab_View<T,D>::const_subview_type
const_Hyperslab_View<T,D>::major_slice(size_type index) const
{
    REQUIRE(index < d_indexer.dims().front());

    // Create slice dimensions
    const auto &subindex = d_indexer.major_slice_dims();

    // Slice the data
    size_type width = std::accumulate(subindex.begin(), subindex.end(), 1u,
                                      std::multiplies<size_type>());
    const_view subview(begin().ptr() + width * index,
                       begin().ptr() + width * (index + 1));

    return const_subview_type(subview, subindex);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access a slice of the hyperslab (most minor index).
 *
 * \warning This returns a strided view field; it will not be contiguous in
 * memory.
 *
 * \warning There will be several integer multiplies and some on-stack integer
 * copies associated with this operation. Save and reuse these slices whenever
 * possible.
 */
template<class T, size_t D>
typename const_Hyperslab_View<T,D>::const_subview_type
const_Hyperslab_View<T,D>::minor_slice(size_type index) const
{
    REQUIRE(index < d_indexer.dims().back());

    // Create slice dimensions
    const auto &subindex = d_indexer.minor_slice_dims();

    // Slice the data
    const_view subview = strided_slice(index,
                                       index + size(),
                                       d_indexer.dims().back());

    return const_subview_type(subview, subindex);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator
 */
template<class T, size_t D>
const_Hyperslab_View<T,D>&
const_Hyperslab_View<T,D>::operator=(const const_Hyperslab_View<T,D>& rhs)
{
    Base::operator=(rhs);
    d_indexer = rhs.d_indexer;
    ENSURE(size() == d_indexer.size());
    return *this;
}

//---------------------------------------------------------------------------//
// LENGTH-ONE SPECIALIZATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Construct with a 1D view.
 */
template<class T>
Hyperslab_View<T,1u>::Hyperslab_View(
        view       data,
        index_type dims)
    : Base(data)
{
    REQUIRE(data.size() == dims.front());
}

//---------------------------------------------------------------------------//
// LENGTH-ONE CONST SPECIALIZATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Construct with a 1D view.
 */
template<class T>
const_Hyperslab_View<T,1u>::const_Hyperslab_View(
        const_view data,
        index_type dims)
    : Base(data)
{
    REQUIRE(data.size() == dims.front());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct with a 1D view.
 */
template<class T>
const_Hyperslab_View<T,1u>::const_Hyperslab_View(
        const Hyperslab_View<T,1u> &view)
    : Base(view)
{
    /* * */
}

} // end namespace profugus

#endif // Utils_utils_Hyperslab_View_i_hh

//---------------------------------------------------------------------------//
//                 end of Hyperslab_View.i.hh
//---------------------------------------------------------------------------//
