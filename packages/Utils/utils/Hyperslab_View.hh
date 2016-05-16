//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/Hyperslab_View.hh
 * \author Seth R Johnson
 * \date   Mon Dec 08 15:25:35 2014
 * \brief  Hyperslab_View class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Hyperslab_View_hh
#define Utils_utils_Hyperslab_View_hh

#include "View_Field.hh"

#include <cstddef>
#include "Regular_Indexer.hh"

namespace profugus
{

// Forward-declare const classes so that they can reference each other
template<class T, std::size_t D> class Hyperslab_View;
template<class T, std::size_t D> class const_Hyperslab_View;

//===========================================================================//
/*!
 * \class Hyperslab_View
 * \brief View field with indexing into multiple dimensions.
 *
 * This is used for multi-dimensional data indexing, allowing inline operations
 * of the form: \code
   Hyperslab_View<double, 3> a(data, dims);

   a[2][1][4] = 2.0;
 \endcode
 *
 * The hyperslab still supports standard field-like operations for begin, end,
 * etc., so that it can be filled, transmitted, and compared.
 *
 * \warning There is a necessary inconsistency between the result of operator[]
 * and the value type. The value and iterators are field-like, accessing single
 * scalar values; but the bracket operator returns
 */
/*!
 * \example utils/test/tstHyperslab_View.cc
 *
 * Test of Hyperslab_View.
 */
//===========================================================================//

template<class T, std::size_t D>
class Hyperslab_View : public View_Field<T>
{
    typedef View_Field<T> Base;
  public:
    //@{
    //! Typedefs
    typedef std::size_t                    size_type;
    typedef Regular_Indexer<size_type, D>  Indexer_t;
    typedef typename Indexer_t::index_type index_type;

    typedef View_Field<T>            view;
    typedef typename view::reference reference;
    typedef Hyperslab_View<T, D-1>   subview_type;

    typedef const_View_Field<T>                  const_view;
    typedef typename const_view::const_reference const_reference;
    typedef const_Hyperslab_View<T, D-1>         const_subview_type;
    //@}

  public:
    // >>> CONSTRUCTOR

    // Construct with empty dimensions
    inline Hyperslab_View();

    // Construct with data pointers, dimensionality
    inline Hyperslab_View(view data, index_type dims);

    // >>> ACCESSORS

    // Access a single value
    inline reference operator[] (const index_type& index);
    inline const_reference operator[] (const index_type& index) const;

    // Slice the most major axis
    inline subview_type operator[] (size_type i) { return major_slice(i); }
    inline const_subview_type operator[] (size_type i) const
    { return major_slice(i); }

    // Slice the most major axis
    inline subview_type major_slice(size_type index);
    inline const_subview_type major_slice(size_type index) const;

    // Slice the most minor axis
    inline subview_type minor_slice(size_type index);
    inline const_subview_type minor_slice(size_type index) const;

    //! Access dimensions
    const index_type& dims() const { return d_indexer.dims(); }

    // View field methods
    using Base::front;
    using Base::back;
    using Base::size;
    using Base::stride;
    using Base::empty;

  public:
    // View field methods
    using Base::data;
    using Base::begin;
    using Base::end;
    using Base::cbegin;
    using Base::cend;

    // Assignment
    inline Hyperslab_View<T,D>& operator=(const Hyperslab_View<T,D>& rhs);

    // >>> FRIENDS
    friend class const_Hyperslab_View<T,D>;

  protected:
    // >>> IMPLEMENTATION

    Indexer_t d_indexer;

    using Base::strided_slice;
};

//===========================================================================//
/*!
 * \class const_Hyperslab_View
 * \brief View field with indexing into multiple dimensions.
 *
 * This is used for multi-dimensional data indexing.
 */
//===========================================================================//

template<class T, std::size_t D>
class const_Hyperslab_View : public const_View_Field<T>
{
    typedef const_View_Field<T> Base;
  public:
    //@{
    //! Typedefs
    typedef std::size_t                     size_type;
    typedef Regular_Indexer<size_type, D> Indexer_t;
    typedef typename Indexer_t::index_type  index_type;

    typedef const_View_Field<T>                  const_view;
    typedef typename const_view::const_reference const_reference;
    typedef const_Hyperslab_View<T, D-1>         const_subview_type;
    //@}

  public:
    // >>> CONSTRUCTOR

    // Construct with empty view, no dimensions
    inline const_Hyperslab_View();

    // Construct with data pointers, dimensionality
    inline const_Hyperslab_View(const_view data, index_type dims);

    // Construct implicitly from mutable hyperslab
    inline const_Hyperslab_View(const Hyperslab_View<T,D> &view);

    // >>> ACCESSORS

    // Access a single value
    inline const_reference operator[] (const index_type& index) const;

    // Slice the most major axis
    inline const_subview_type operator[] (size_type i) const
    { return major_slice(i); }

    // Slice the most major axis
    inline const_subview_type major_slice(size_type index) const;

    // Slice the most minor axis
    inline const_subview_type minor_slice(size_type index) const;

    //! Access dimensions
    const index_type& dims() const { return d_indexer.dims(); }

    // View field methods
    using Base::front;
    using Base::back;
    using Base::size;
    using Base::stride;
    using Base::empty;

  public:
    // View field methods
    using Base::data;
    using Base::begin;
    using Base::end;
    using Base::cbegin;
    using Base::cend;

    // Assignment
    inline const_Hyperslab_View<T,D>& operator=(
            const const_Hyperslab_View<T,D>& rhs);

  protected:
    // >>> IMPLEMENTATION

    Indexer_t d_indexer;

    using Base::strided_slice;
};

//===========================================================================//
/*!
 * One-element hyperslab is a simple view field.
 *
 * It does not have or need an indexer.
 */
template<class T>
class Hyperslab_View<T, 1u> : public View_Field<T>
{
    typedef View_Field<T> Base;
  public:
    typedef View_Field<T>             view;
    typedef typename view::size_type  size_type;
    typedef Vector_Lite<size_type, 1> index_type;

  public:
    // Construct with the 1-D view
    inline Hyperslab_View(view data, index_type dims);

    // View field methods
    using Base::front;
    using Base::back;
    using Base::size;
    using Base::stride;
    using Base::empty;

  public:
    // View field methods
    using Base::data;
    using Base::begin;
    using Base::end;
    using Base::cbegin;
    using Base::cend;
};

//===========================================================================//
/*!
 * One-element hyperslab is a simple view field.
 *
 * It does not have or need an indexer.
 */
template<class T>
class const_Hyperslab_View<T, 1u> : public const_View_Field<T>
{
    typedef const_View_Field<T> Base;
  public:
    typedef const_View_Field<T>            const_view;
    typedef typename const_view::size_type size_type;
    typedef Vector_Lite<size_type, 1>      index_type;

  public:
    // Construct with the 1-D view
    inline const_Hyperslab_View(const_view data, index_type dims);

    // Construct implicitly from a non-const view field
    inline const_Hyperslab_View(const Hyperslab_View<T,1u> &view);

    // View field methods
    using Base::front;
    using Base::back;
    using Base::size;
    using Base::stride;
    using Base::empty;

  public:
    // View field methods
    using Base::data;
    using Base::begin;
    using Base::end;
    using Base::cbegin;
    using Base::cend;
};

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE METHODS
//---------------------------------------------------------------------------//
#include "Hyperslab_View.i.hh"
//---------------------------------------------------------------------------//

#endif // Utils_utils_Hyperslab_View_hh

//---------------------------------------------------------------------------//
//                 end of Hyperslab_View.hh
//---------------------------------------------------------------------------//
