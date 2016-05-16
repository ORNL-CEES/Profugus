//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/View_Field_Iterator.i.hh
 * \author Gregory G. Davidson
 * \date   Wed Oct 08 11:07:08 2014
 * \brief  Member definitions of class View_Field_Iterator.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_View_Field_Iterator_i_hh
#define Utils_utils_View_Field_Iterator_i_hh

namespace profugus
{

//===========================================================================//
// VF_ITERATOR IMPLEMENTATION
//===========================================================================//

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
template<typename T>
VF_Iterator<T>::VF_Iterator()
    : d_ptr(0)
    , d_stride(1)
{   }

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename T>
VF_Iterator<T>::VF_Iterator(pointer      ptr_in,
                            unsigned int stride)
    : d_ptr(ptr_in)
    , d_stride(stride)
{
    REQUIRE(stride > 0);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Equality operator
 */
template<typename T>
bool VF_Iterator<T>::operator==(const This &iter) const
{
    REQUIRE(stride() == iter.stride());
    return d_ptr == iter.d_ptr;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Inequality operator
 */
template<typename T>
bool VF_Iterator<T>::operator!=(const This &iter) const
{
    REQUIRE(stride() == iter.stride());
    return d_ptr != iter.d_ptr;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Postfix increment.
 */
template<typename T>
VF_Iterator<T> VF_Iterator<T>::operator++(int)
{
    This tmp = *this;
    d_ptr += d_stride;
    return tmp;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Postfix decrement.
 */
template<typename T>
VF_Iterator<T> VF_Iterator<T>::operator--(int)
{
    This tmp = *this;
    d_ptr -= d_stride;
    return tmp;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Addition operator
 */
template<typename T>
VF_Iterator<T> VF_Iterator<T>::operator+(difference_type n) const
{
    return This(d_ptr + n * d_stride, d_stride);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compound addition
 */
template<typename T>
VF_Iterator<T>& VF_Iterator<T>::operator+=(difference_type n)
{
    d_ptr += n * d_stride;
    return *this;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Subtraction operator.
 */
template<typename T>
VF_Iterator<T> VF_Iterator<T>::operator-(difference_type n) const
{
    return This(d_ptr - n * d_stride, d_stride);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Subtraction operator.
 */
template<typename T>
typename VF_Iterator<T>::difference_type
VF_Iterator<T>::operator-(const This &vf_iter) const
{
    difference_type delta = d_ptr - vf_iter.d_ptr;
    if (d_stride == 1)
        return delta;
    return delta / d_stride;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Offset dereference operator.
 */
template<typename T>
typename VF_Iterator<T>::reference
VF_Iterator<T>::operator[](int n)
{
    REQUIRE(d_ptr);

    return d_ptr[n*d_stride];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constant offset dereference operator.
 */
template<typename T>
const typename VF_Iterator<T>::reference
VF_Iterator<T>::operator[](int n) const
{
    REQUIRE(d_ptr);

    return d_ptr[n*d_stride];
}

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Addition operator.
 */
template<typename T>
VF_Iterator<T> operator+(
        typename VF_Iterator<T>::difference_type  n,
        const VF_Iterator<T>                     &iter)
{
    return VF_Iterator<T>(iter.get_pointer + n*iter.stride(), iter.stride());
}

//===========================================================================//
// CONST_VF_ITERATOR IMPLEMENTATION
//===========================================================================//

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
template<typename T>
const_VF_Iterator<T>::const_VF_Iterator()
    : d_ptr(0)
    , d_stride(1)
{   }

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename T>
const_VF_Iterator<T>::const_VF_Iterator(pointer ptr_in, unsigned int stride)
    : d_ptr(ptr_in)
    , d_stride(stride)
{
    REQUIRE(stride > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor from a VF_Iterator.
 */
template<typename T>
const_VF_Iterator<T>::const_VF_Iterator(const VF_Iterator<T> &iter)
    : d_ptr(iter.ptr())
    , d_stride(iter.stride())
{
    ENSURE(d_stride > 0);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Equality operator
 */
template<typename T>
bool const_VF_Iterator<T>::operator==(const This &iter) const
{
    REQUIRE(stride() == iter.stride());
    return d_ptr == iter.d_ptr;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Inequality operator
 */
template<typename T>
bool const_VF_Iterator<T>::operator!=(const This &iter) const
{
    REQUIRE(stride() == iter.stride());
    return d_ptr != iter.d_ptr;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Postfix increment.
 */
template<typename T>
const_VF_Iterator<T> const_VF_Iterator<T>::operator++(int)
{
    This tmp = *this;
    d_ptr += d_stride;
    return tmp;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Postfix decrement.
 */
template<typename T>
const_VF_Iterator<T> const_VF_Iterator<T>::operator--(int)
{
    This tmp = *this;
    d_ptr -= d_stride;
    return tmp;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Addition operator
 */
template<typename T>
const_VF_Iterator<T> const_VF_Iterator<T>::operator+(difference_type n) const
{
    return This(d_ptr + n * d_stride, d_stride);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compound addition
 */
template<typename T>
const_VF_Iterator<T>& const_VF_Iterator<T>::operator+=(difference_type n)
{
    d_ptr += n * d_stride;
    return *this;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Subtraction operator.
 */
template<typename T>
const_VF_Iterator<T> const_VF_Iterator<T>::operator-(difference_type n) const
{
    return This(d_ptr - n * d_stride, d_stride);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Subtraction operator.
 */
template<typename T>
typename const_VF_Iterator<T>::difference_type
const_VF_Iterator<T>::operator-(const This &vf_iter) const
{
    difference_type delta = d_ptr - vf_iter.d_ptr;
    if (d_stride == 1)
        return delta;
    return delta / d_stride;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Offset dereference operator.
 */
template<typename T>
typename const_VF_Iterator<T>::reference
const_VF_Iterator<T>::operator[](int n) const
{
    REQUIRE(d_ptr != 0);

    return d_ptr[n*d_stride];
}

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Addition operator.
 */
template<typename T>
const_VF_Iterator<T> operator+(
        typename const_VF_Iterator<T>::difference_type  n,
        const const_VF_Iterator<T>                     &iter)
{
    return const_VF_Iterator<T>(iter.get_pointer + n*iter.stride(),
                                iter.stride());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Inequality operator between const/non-const VF
 */
template<typename T>
bool operator==(const VF_Iterator<T>       &left,
                const const_VF_Iterator<T> &right)
{
    return left.ptr() == right.ptr();
}

} // end namespace profugus

#endif // Utils_utils_View_Field_Iterator_i_hh

//---------------------------------------------------------------------------//
//                 end of View_Field_Iterator.i.hh
//---------------------------------------------------------------------------//
