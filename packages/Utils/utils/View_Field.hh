//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Nemesis/utils/View_Field.hh
 * \author Gregory G. Davidson
 * \date   Thu Aug  9 09:47:05 2007
 * \brief  View_Field class declaration.
 * \note   Copyright (C) 2007 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Nemesis_utils_View_Field_hh
#define Nemesis_utils_View_Field_hh

#include <cstring>
#include <cstddef>
#include <iterator>
#include <vector>

#include "Nemesis/harness/DBC.hh"
#include "Nemesis/harness/SWIG.hh"
#include "View_Field_Iterator.hh"

namespace nemesis
{

//===========================================================================//
/*!
 * \class View_Field
 * \brief Provides views into field data.
 *
 * This class provides a common wrapper around random access containers,
 * and is able to return slices.
 *
 * \note  The underlying data \b must be contiguous.
 */
//===========================================================================//

template<typename T>
class View_Field
{
  public:
    //@{
    //! Container typedefs
    typedef T        value_type;
    typedef T&       reference;
    typedef const T& const_reference;
    typedef T*       pointer;
    typedef const T* const_pointer;

    typedef VF_Iterator<value_type>       iterator;
    typedef const_VF_Iterator<value_type> const_iterator;

    typedef std::reverse_iterator<iterator>       reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    typedef std::size_t    size_type;
    typedef std::ptrdiff_t difference_type;
    typedef unsigned int   stride_type;

    typedef View_Field<T> View_Field_t;
    //@}

  private:
    // >>> PRIVATE DATA
    // Pointer to the beginning of the field data
    iterator d_begin_iterator;
    // Pointer to the ending of the field data
    iterator d_end_iterator;

  SWIG_PRIVATE:
    // Default constructor.  Sets both iterators at null
    inline View_Field();

    // Constructor. Takes two pointers and an optional stride parameter
    inline View_Field(pointer     first,
                      pointer     last,
                      stride_type stride = 1);

    // Construct from view field iterators
    inline View_Field(iterator     first,
                      iterator     last);

    // Assign operator
    // This should be used in favor of std::copy, since it has an optimization
    // for stride 1 view fields
    void assign(const_iterator first,
                const_iterator last);
    // General assign operator
    template<typename InputIterator>
    void assign(InputIterator first,
                InputIterator last);

  public:
    // ACCESSORS

    // Returns field value at index i.
    inline reference operator[](size_type i);

    // Returns constant field value at index i.
    inline const_reference operator[](size_type i) const;

    // Return the front piece of data
    inline reference front();
    inline const_reference front() const;

    // Return the back piece of data
    inline reference back();
    inline const_reference back() const;

    // Number of field values.
    inline size_type size() const;

    // Stride between elements
    inline stride_type stride() const;

    //! Returns true if the size of the field is zero.
    bool empty() const { return d_begin_iterator == d_end_iterator; }

  SWIG_PRIVATE:

    // >>> RAW DATA ACCESS (like C++11)
    //@{
    //! Raw data access
    inline pointer data();
    inline const_pointer data() const;
    //@}

    // >>> ITERATORS

    //! Begin iterator.
    iterator begin() { return d_begin_iterator; }

    //@{
    //! Begin const_iterator.
    const_iterator begin() const  { return d_begin_iterator; }
    const_iterator cbegin() const { return d_begin_iterator; }
    //@}

    //! End iterator.
    iterator end() { return d_end_iterator; }

    //@{
    //! End const_iterator.
    const_iterator end() const  { return d_end_iterator; }
    const_iterator cend() const { return d_end_iterator; }
    //@}

    // >>> REVERSE ITERATORS

    //! Begin reverse iterator.
    reverse_iterator rbegin() { return reverse_iterator(end()); }

    //@{
    //! Begin reverse const_iterator.
    inline const_reverse_iterator rbegin() const;
    const_reverse_iterator crbegin() const { return rbegin(); }
    //@}

    //! End reverse iterator.
    reverse_iterator rend() { return reverse_iterator(begin()); }

    //@{
    //! End reverse const_iterator.
    inline const_reverse_iterator rend() const;
    const_reverse_iterator crend() const { return rend(); }
    //@}

  public:
    // Divide the field into slices of size slice_size, and return the
    // slice_num-th slice.
    inline View_Field_t slice(size_type slice_size,
                              int       slice_num) const;

    // Divide the field into slices of size slice_size, and return slices
    // slice_begin to slice_end.
    inline View_Field_t slice(size_type slice_size,
                              int       slice_begin,
                              int       slice_end) const;

    // Return a slice from [index_begin, index_end).
    inline View_Field_t general_slice(size_type index_begin,
                                      size_type index_end) const;

    // Return a strided slice from [index_begin, index_end)
    inline View_Field_t strided_slice(size_type index_begin,
                                      size_type index_end,
                                      size_type stride) const;
  private:
    // For bounds checking.
    inline bool valid_index(const size_type i) const;
};

//===========================================================================//
/*!
 * \class const_View_Field
 *
 * \brief Provides constant views into the field data.
 *
 * This class provides a common wrapper around random access containers,
 * especially State fields and Epetra_Vectors holding flux moments.  With a
 * const_View_Field, the data may not be changed.
 */
//===========================================================================//

template<typename T>
class const_View_Field
{
  public:
    //@{
    //! Typedefs
    typedef T        value_type;
    typedef const T& const_reference;
    typedef const T* const_pointer;

    typedef const_VF_Iterator<value_type>         const_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    typedef std::size_t    size_type;
    typedef std::ptrdiff_t difference_type;
    typedef unsigned int   stride_type;

    typedef View_Field<T>       View_Field_t;
    typedef const_View_Field<T> const_View_Field_t;
    //@}

  private:
    /// >>> PRIVATE DATA
    const_iterator d_begin_iterator;
    const_iterator d_end_iterator;

  SWIG_PRIVATE:
    // Default constructor.
    inline const_View_Field();

    // Constructor.  Build implicitly from a View_Field.
    inline const_View_Field(const View_Field_t& field_view_in);

    // Constructor. Takes two pointers and an optional stride
    inline const_View_Field(const_pointer first,
                            const_pointer last,
                            stride_type   stride = 1);

    // Construct from view field iterators
    inline const_View_Field(const_iterator     first,
                            const_iterator     last);

  public:
    // ACCESSORS

    // Returns constant field value at index i.
    inline const_reference operator[](const size_type i) const;

    // Return the front data
    inline const_reference front() const;

    // Return the back data
    inline const_reference back() const;

    // Number of field values.
    inline size_type size() const;

    // Stride between elements
    inline stride_type stride() const;

    //! Returns true if the size of the field is zero.
    bool empty() const { return d_begin_iterator == d_end_iterator; }

  SWIG_PRIVATE:
    // >>> RAW DATA ACCESS (like C++11)
    //! Raw data access
    inline const_pointer data() const;

    // >>> ITERATORS
    //@{
    //! Begin const_iterator.
    const_iterator begin() const  { return d_begin_iterator; }
    const_iterator cbegin() const { return d_begin_iterator; }
    //@}

    //@{
    //! End const_iterator.
    const_iterator end() const  { return d_end_iterator; }
    const_iterator cend() const { return d_end_iterator; }
    //@}

    // >>> REVERSE ITERATORS

    //@{
    //! Begin reverse const_iterator.
    inline const_reverse_iterator rbegin() const;
    const_reverse_iterator crbegin() const { return rbegin(); }
    //@}

    //@{
    //! End reverse const_iterator.
    inline const_reverse_iterator rend() const;
    const_reverse_iterator crend() const { return rend(); }
    //@}

  public:
    // Divide the view into slices of size slice_size and return the
    // slice_num-th slice.
    inline const_View_Field_t slice(size_type slice_size,
                                    int       slice_num) const;

    // Divide the field into slices of size slice_size, and return slices
    // slice_begin to slice_end.
    inline const_View_Field_t slice(size_type slice_size,
                                    int       slice_begin,
                                    int       slice_end) const;

    // Return a slice from [index_begin, index_end).
    inline const_View_Field_t general_slice(size_type index_begin,
                                            size_type index_end) const;

    // Return a strided slice from [index_begin, index_end)
    inline const_View_Field_t strided_slice(size_type index_begin,
                                            size_type index_end,
                                            size_type stride) const;

  private:
    // For bounds checking.
    inline bool valid_index(const size_type i) const;
};

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// View of raw arrays

template<typename T, std::size_t N>
inline View_Field<T> make_view(T (&array)[N]);

template<typename T, std::size_t N>
inline const_View_Field<T> make_view(const T (&array)[N]);

template<typename T, std::size_t N>
inline const_View_Field<T> make_const_view(const T (&array)[N]);

// View of vectors

template<typename T, class A>
inline View_Field<T> make_view(std::vector<T,A> &vec);

template<typename T, class A>
inline const_View_Field<T> make_view(const std::vector<T,A> &vec);

template<typename T, class A>
inline const_View_Field<T> make_const_view(const std::vector<T,A> &vec);

// Views of views

template<typename T>
inline View_Field<T> make_view(nemesis::View_Field<T> view);

template<typename T>
inline const_View_Field<T> make_view(nemesis::const_View_Field<T> view);

template<typename T>
inline const_View_Field<T> make_const_view(nemesis::const_View_Field<T> view);


//---------------------------------------------------------------------------//

} // end namespace nemesis

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "View_Field.i.hh"

#endif // Nemesis_utils_View_Field_hh

//---------------------------------------------------------------------------//
//              end of nemesis/View_Field.hh
//---------------------------------------------------------------------------//
