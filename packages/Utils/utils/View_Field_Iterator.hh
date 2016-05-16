//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/View_Field_Iterator.hh
 * \author Gregory G. Davidson
 * \date   Wed Oct 08 11:07:08 2014
 * \brief  View_Field_Iterator class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_View_Field_Iterator_hh
#define Utils_utils_View_Field_Iterator_hh

#include <iterator>

#include "Utils/harness/DBC.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class VF_Iterator
 * \brief Provides a mutable iterator with an optional stride.
 */
/*!
 * \example utils/test/tstView_Field_Iterator.cc
 *
 * Test of VF_Iterator.
 */
//===========================================================================//

template<typename T>
class VF_Iterator : public std::iterator_traits<T*>
{
  public:
    //@{
    //! Useful typedefs
    typedef std::iterator_traits<T*>       Base;
    typedef typename Base::pointer         pointer;
    typedef typename Base::value_type      value_type;
    typedef unsigned int                   stride_type;
    typedef typename Base::reference       reference;
    typedef typename Base::difference_type difference_type;
    typedef VF_Iterator<T>                 This;
    //@}

  private:
    // Data
    pointer d_ptr;
    stride_type d_stride;

  public:
    // Default constructor
    inline VF_Iterator();

    // Constructor
    inline explicit VF_Iterator(pointer ptr_in, stride_type stride);

    //! Return the stride
    stride_type stride() const { return d_stride; }
    //! Return the underlying pointer
    pointer ptr() const { return d_ptr; }
    //! Return whether the pointer points to valid memory
    bool is_valid() const { return d_ptr != 0; }

    //@{
    //! Dereference operator.
    reference operator*()             { REQUIRE(d_ptr != 0); return *d_ptr; }
    const reference operator*() const { REQUIRE(d_ptr != 0); return *d_ptr; }
    //@}

    //@{
    //! Pointer operator.
    pointer operator->()             { REQUIRE(d_ptr != 0); return d_ptr; }
    const pointer operator->() const { REQUIRE(d_ptr != 0); return d_ptr; }
    //@}

    // Equality operator with pointer
    inline bool operator==(const This& iter) const;
    // Inequality operator with pointer
    inline bool operator!=(const This& iter) const;
    // Non-null operator
#ifdef __NVCC__
    operator void*() const { return d_ptr; }
#else
    explicit operator bool() const { return d_ptr != nullptr; }
#endif

    //! Prefix increment
    This& operator++() { d_ptr += d_stride; return *this; }
    // Postfix increment
    inline This operator++(int);

    //! Prefix decrement
    This& operator--() { d_ptr -= d_stride; return *this; }
    // Postfix decrement
    inline This operator--(int);

    // Addition operator
    inline This operator+(difference_type n) const;
    // Compound addition
    This& operator+=(difference_type n);

    //! Subtraction operator
    inline This operator-(difference_type n) const;
    inline difference_type operator-(const This &vf_iter) const;

    //! Compound subtraction
    This& operator-=(difference_type n) { d_ptr -= n*d_stride; return *this; }

    //@{
    //! Comparators.
    bool operator<(const This &iter) const  { return d_ptr < iter.d_ptr; }
    bool operator>(const This &iter) const  { return d_ptr > iter.d_ptr; }
    bool operator<=(const This &iter) const { return d_ptr <= iter.d_ptr; }
    bool operator>=(const This &iter) const { return d_ptr >= iter.d_ptr; }
    //@}

    // Offset dereference operator.
    inline reference operator[](int n);
    inline const reference operator[](int n) const;
};

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Addition operator.
 */
template<typename T>
inline VF_Iterator<T> operator+(
        typename VF_Iterator<T>::difference_type  n,
        const VF_Iterator<T>                     &iter);

//===========================================================================//
/*!
 * \class const_VF_Iterator
 * \brief Provides a constant iterator with an optional stride.
 */
/*!
 * \example utils/test/tstView_Field_Iterator.cc
 *
 * Test of const_VF_Iterator.
 */
//===========================================================================//

template<typename T>
class const_VF_Iterator : public std::iterator_traits<const T*>
{
  public:
    //@{
    //! Useful typedefs
    typedef std::iterator_traits<const T*> Base;
    typedef typename Base::pointer         pointer;
    typedef typename Base::value_type      value_type;
    typedef unsigned int                   stride_type;
    typedef typename Base::difference_type difference_type;
    typedef typename Base::reference       reference;
    typedef const_VF_Iterator<T>           This;
    //@}

  private:
    // Data
    pointer d_ptr;
    stride_type d_stride;

  public:
    // Default constructor
    inline const_VF_Iterator();

    // Constructor
    inline const_VF_Iterator(const pointer ptr_in, stride_type stride);

    // Constructor from a VF_Iterator
    inline const_VF_Iterator(const VF_Iterator<T> &iter);

    //! Return the stride.
    stride_type stride() const { return d_stride; }
    //! Return the underlying pointer.
    pointer ptr() const { return d_ptr; }
    //! Return whether the pointer points to valid memory.
    bool is_valid() const { return d_ptr != 0; }

    //! Dereference operator.
    reference operator*() const { REQUIRE(d_ptr != 0); return *d_ptr; }

    //! Pointer operator.
    pointer operator->() const { REQUIRE(d_ptr != 0); return d_ptr; }

    // Equality operator
    inline bool operator==(const This &iter) const;
    // Inequality operator
    inline bool operator!=(const This &iter) const;
    // Non-null operator
#ifdef __NVCC__
    operator void*() const { return d_ptr; }
#else
    explicit operator bool() const { return d_ptr != nullptr; }
#endif

    //! Prefix increment
    This& operator++() { d_ptr += d_stride; return *this; }
    // Postfix increment
    inline This operator++(int);

    //! Prefix decrement
    This& operator--() { d_ptr -= d_stride; return *this; }
    // Postfix decrement
    inline This operator--(int);

    // Addition operator
    inline This operator+(difference_type n) const;
    // Compound addition
    This& operator+=(difference_type n);

    //! Subtraction operator
    inline This operator-(difference_type n) const;
    inline difference_type operator-(const This &vf_iter) const;

    //! Compound subtraction
    This& operator-=(difference_type n) { d_ptr -= n*d_stride; return *this; }

    //@{
    //! Comparators.
    bool operator<(const This &iter) const  { return d_ptr < iter.d_ptr; }
    bool operator>(const This &iter) const  { return d_ptr > iter.d_ptr; }
    bool operator<=(const This &iter) const { return d_ptr <= iter.d_ptr; }
    bool operator>=(const This &iter) const { return d_ptr >= iter.d_ptr; }
    //@}

    // Offset dereference operator.
    inline reference operator[](int n) const;
};

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Addition operator.
 */
template<typename T>
inline const_VF_Iterator<T> operator+(
        typename const_VF_Iterator<T>::difference_type  n,
        const const_VF_Iterator<T>                     &iter);

//---------------------------------------------------------------------------//
/*!
 * \brief Equality operator between const/non-const VF
 */
template<typename T>
inline bool operator==(const VF_Iterator<T>       &left,
                       const const_VF_Iterator<T> &right);

template<typename T>
bool operator==(const const_VF_Iterator<T> &left,
                const VF_Iterator<T>       &right)
{
    return (right == left);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Inequality operator between const/non-const VF
 */
template<typename T>
bool operator!=(const VF_Iterator<T> &left,
                const const_VF_Iterator<T> &right)
{
    return !(left == right);
}

template<typename T>
bool operator!=(const const_VF_Iterator<T> &left,
                const VF_Iterator<T>      &right)
{
    return !(right == left);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE FUNCTION DEFINITIONS
//---------------------------------------------------------------------------//

#include "View_Field_Iterator.i.hh"

#endif // Utils_utils_View_Field_Iterator_hh

//---------------------------------------------------------------------------//
//                 end of View_Field_Iterator.hh
//---------------------------------------------------------------------------//
