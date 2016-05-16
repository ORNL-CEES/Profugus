//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/utils/Hyperslab_Vector.i.hh
 * \author Seth R Johnson
 * \date   Mon Aug 10 15:59:32 2015
 * \brief  Hyperslab_Vector inline method definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Hyperslab_Vector_i_hh
#define Utils_utils_Hyperslab_Vector_i_hh

namespace profugus
{
//---------------------------------------------------------------------------//
/*!
 * \brief Construct an empty vector
 */
template<class T, size_t D>
Hyperslab_Vector<T,D>::Hyperslab_Vector()
{
    ENSURE(empty());
    ENSURE(d_data.size() == d_indexer.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy from a view to other data
 */
template<class T, size_t D>
Hyperslab_Vector<T,D>::Hyperslab_Vector(index_type dims,
                                        value_type val)
    : d_indexer(dims)
    , d_data(d_indexer.size(), val)
{
    ENSURE(d_indexer.dims() == dims);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy from a view to other data
 */
template<class T, size_t D>
Hyperslab_Vector<T,D>::Hyperslab_Vector(const_View_t view)
    : d_indexer(view.dims())
    , d_data(view.begin(), view.end())
{
    ENSURE(d_data.size() == view.size());
    ENSURE(d_data.size() == d_indexer.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Resize to the zero vector
 */
template<class T, size_t D>
void Hyperslab_Vector<T,D>::clear()
{
    d_indexer = Indexer_t();
    d_data.clear();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Resize to a new set of dimensions
 *
 * \warning If you resize any but the longest-stride dimension, existing values
 * in the array will be shifted unexpectedly!
 */
template<class T, size_t D>
void Hyperslab_Vector<T,D>::resize(index_type dims, value_type val)
{
    d_indexer = Indexer_t(dims);
    d_data.resize(d_indexer.size(), val);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign all-new values
 */
template<class T, size_t D>
void Hyperslab_Vector<T,D>::assign(index_type dims, value_type val)
{
    d_indexer = Indexer_t(dims);
    d_data.assign(d_indexer.size(), val);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Swap with another instance
 */
template<class T, size_t D>
void Hyperslab_Vector<T,D>::swap(This& rhs)
{
    using std::swap;
    swap(d_indexer, rhs.d_indexer);
    swap(d_data, rhs.d_data);
}

//---------------------------------------------------------------------------//
} // end namespace profugus

#endif // Utils_utils_Hyperslab_Vector_i_hh

//---------------------------------------------------------------------------//
// end of Utils/utils/Hyperslab_Vector.i.hh
//---------------------------------------------------------------------------//
