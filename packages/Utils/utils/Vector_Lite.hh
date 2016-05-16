//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/Vector_Lite.hh
 * \author Derived from Rob Lowrie's Vector_Lite (LANL)
 * \date   Thu Jan  3 11:39:29 2008
 * \brief  Vector_Lite class definition.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Vector_Lite_hh
#define Utils_utils_Vector_Lite_hh

#include <iostream>
#include <numeric>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <initializer_list>

#include "Utils/harness/DBC.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Vector_Lite
 * \brief Array container that is a wrapper around a standard C array.
 *
 * It adds iterator and arithemtic support, along with bounds checking (via
 * Utils' DBC).
 *
 * An alternative to this class is boost::array (www.boost.org).  However,
 * boost::array is an aggregate type, which has advantages (can use
 * initializers) and disadvantages (public data, cannot be a base class).
 * boost::array also doesn't do bounds checking.
 *
 * \param T Type of each array element.
 *
 * \param N Length of array.
 */
/*!
 * \example utils/test/tstVector_Lite.cc
 *
 * Test of Vector_Lite.
 */
//===========================================================================//

template <class T, size_t N>
class Vector_Lite
{
  public:
    //@{
    //! Typedefs.
    typedef T                 value_type;
    typedef value_type*       pointer;
    typedef const T*          const_pointer;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef ptrdiff_t         difference_type;
    typedef size_t            size_type;
    typedef pointer           iterator;
    typedef const_pointer     const_iterator;

    typedef std::reverse_iterator<iterator>       reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    //@}

  private:
    // >>> DATA

    T d_U[N];

  public:
    // Empty constructor
    inline explicit Vector_Lite();

    // Constructor based on a scalar value.
    inline explicit Vector_Lite(const T &u);

    // Constructor for N = 2.
    inline Vector_Lite(const T &u0, const T &u1);

    // Constructor for N = 3.
    inline Vector_Lite(const T &u0, const T &u1, const T &u2);

    // Constructor for N = 4.
    inline Vector_Lite(const T &u0, const T &u1, const T &u2, const T &u3);

    // Constructor for N = 5.
    inline Vector_Lite(const T &u0, const T &u1, const T &u2,
                       const T &u3, const T &u4);

    // Initializer list construction.
    inline Vector_Lite(std::initializer_list<T> list);

    // Copy and move constructor
    Vector_Lite(const Vector_Lite& rhs) = default;
    Vector_Lite(Vector_Lite&& rhs) = default;

    // >>> MANIPULATORS

    // Assignment from a another Vector_Lite.
    Vector_Lite &operator=(const Vector_Lite &rhs) = default;
    Vector_Lite &operator=(Vector_Lite &&rhs) = default;

    // Assignment from a scalar.
    inline Vector_Lite &operator=(const T &rhs);

    // Comparisons to another Vector_Lite.
    inline bool operator==(const Vector_Lite &a) const;
    inline bool operator<(const Vector_Lite &a) const;
    inline bool all_gt(const Vector_Lite &a) const;
    inline bool all_lt(const Vector_Lite &a) const;
    inline bool all_ge(const Vector_Lite &a) const;
    inline bool all_le(const Vector_Lite &a) const;

    // Basic arithmetic operations, vector right-hand side.
    inline Vector_Lite &operator+=(const Vector_Lite &a);
    inline Vector_Lite &operator-=(const Vector_Lite &a);
    inline Vector_Lite &operator*=(const Vector_Lite &a);
    inline Vector_Lite &operator/=(const Vector_Lite &a);

    // Basic arithmetic operations, scalar right-hand side.
    template<class T2>
    inline Vector_Lite &operator+=(const T2 &a);
    template<class T2>
    inline Vector_Lite &operator-=(const T2 &a);
    template<class T2>
    inline Vector_Lite &operator*=(const T2 &a);
    template<class T2>
    inline Vector_Lite &operator/=(const T2 &a);

    // >>> ACCESSORS

    //! Indexing using [].
    reference operator[](const size_type i)
    {
        REQUIRE(i < N);
        return d_U[i];
    }

    //! const indexing using [].
    const_reference operator[](const size_type i) const
    {
        REQUIRE(i < N);
        return d_U[i];
    }

    //@{
    //! C++11-like data access
    pointer data() { return d_U; }
    const_pointer data() const { return d_U; }
    //@}

    //! Front and back elements (unchecked)
    const_reference front() const { return d_U[0]; }
    reference front() { return d_U[0]; }
    const_reference back() const { return d_U[N-1]; }
    reference back() { return d_U[N-1]; }
    //@}

    // >>> ITERATOR SUPPORT

    //! Iterator begin.
    iterator begin() { return d_U; }

    //! Const iterator begin.
    const_iterator begin() const { return d_U; }

    //! Begin reverse iterator.
    reverse_iterator rbegin() { return reverse_iterator(end()); }

    //! Begin reverse const_iterator.
    const_reverse_iterator rbegin() const
        { return const_reverse_iterator(end()); }
    const_reverse_iterator crbegin() const { return rbegin(); }

    //! Iterator end.
    iterator end() { return d_U + N; }

    //! End reverse iterator.
    reverse_iterator rend() { return reverse_iterator(begin()); }

    //! End reverse const_iterator.
    const_reverse_iterator rend() const
        { return const_reverse_iterator(begin()); }
    const_reverse_iterator crend() const { return rend(); }

    //! Const iterator end.
    const_iterator end() const { return d_U + N; }

    //! Number of elements (\a N); for STL support.
    size_type size() const { return N; }

    //! Max number of elements (\a N); for STL support.
    size_type max_size() const { return N; }

    //! True if \a N = 0; for STL support.
    bool empty() const { return N == 0; }
};

} // end namespace profugus
//---------------------------------------------------------------------------//
// INLINE AND FREE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Vector_Lite.i.hh"

#endif // Utils_utils_Vector_Lite_hh

//---------------------------------------------------------------------------//
//              end of utils/Vector_Lite.hh
//---------------------------------------------------------------------------//
