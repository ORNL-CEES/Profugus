//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/utils/Hyperslab_Vector.hh
 * \author Seth R Johnson
 * \date   Mon Aug 10 15:59:32 2015
 * \brief  Hyperslab_Vector class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Hyperslab_Vector_hh
#define Utils_utils_Hyperslab_Vector_hh

#include <cstddef>
#include <vector>

#include "Regular_Indexer.hh"
#include "Hyperslab_View.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Hyperslab_Vector
 * \brief Self-contained vector with multiple dimensions.
 *
 * This is sort of like a hyperslab view that contains its own data. It is a
 * simple class meant for local storage. You should access it through the
 * views.
 */
/*!
 * \example utils/test/tstHyperslab_Vector.cc
 *
 * Test of Hyperslab_Vector.
 */
//===========================================================================//

template<class T, std::size_t D>
class Hyperslab_Vector
{
    typedef Hyperslab_Vector<T, D> This;
  public:
    //@{
    //! Typedefs
    typedef Hyperslab_View<T, D>           View_t;
    typedef const_Hyperslab_View<T, D>     const_View_t;
    typedef std::size_t                    size_type;
    typedef Regular_Indexer<size_type,D>   Indexer_t;
    typedef typename Indexer_t::index_type index_type;

    typedef std::vector<T> Storage_t;
    typedef T              value_type;
    //@}

  private:
    // >>> DATA

    // Indexer (multi-dimensions)
    Indexer_t d_indexer;

    // Data storage
    Storage_t d_data;

  public:

    // Construct an empty vector
    inline Hyperslab_Vector();

    // Construct with a given size and fill value
    explicit inline Hyperslab_Vector(index_type dims,
                                     value_type val = value_type());

    // Copy from a view to other data
    explicit inline Hyperslab_Vector(const_View_t view);

    // Resize to the zero vector
    inline void clear();

    // Resize to a new set of dimensions
    inline void resize(index_type dims, value_type val = value_type());

    // Assign a new size and values
    inline void assign(index_type dims, value_type val = value_type());

    // Swap with another one of us (cheap)
    inline void swap(This& rhs);

    // >>> ACCESSORS

    //! Whether the size is zero
    bool empty() const { return d_data.empty(); }

    //! The number of elements
    size_type size() const { return d_data.size(); }

    //! Access indexer
    const Indexer_t& indexer() const { return d_indexer; }

    //! Access dimensions
    const index_type& dims() const { return d_indexer.dims(); }

    //! Access raw data
    //@{
    const value_type* data() const { return d_data.data(); }
    value_type* data() { return d_data.data(); }
    //@}

    //! Access a mutable view
    View_t view()
    {
        return View_t(profugus::make_view(d_data), dims());
    }

    //@{
    //! Access a const view
    const_View_t view() const
    {
        return const_View_t(profugus::make_view(d_data), dims());
    }
    const_View_t cview() const { return view(); }
    //@}

};

//---------------------------------------------------------------------------//
// INLINE TEMPLATE FUNCTIONS
//---------------------------------------------------------------------------//
//! Cheap swap operation
template <class T, std::size_t D>
inline void swap(Hyperslab_Vector<T,D>& a, Hyperslab_Vector<T,D>& b)
{
    a.swap(b);
}

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "Hyperslab_Vector.i.hh"
//---------------------------------------------------------------------------//
#endif // Utils_utils_Hyperslab_Vector_hh

//---------------------------------------------------------------------------//
// end of Utils/utils/Hyperslab_Vector.hh
//---------------------------------------------------------------------------//
