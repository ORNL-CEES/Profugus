//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Vector_Lite.i.hh
 * \author Thomas M. Evans
 * \date   Thu Jan  3 11:39:29 2008
 * \brief  Member definitions of class Vector_Lite.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Vector_Lite_i_hh
#define utils_Vector_Lite_i_hh

namespace profugus
{

//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor based on a scalar value.
 *
 * Initializes all values to \a u.  Note that this ctor also acts as the
 * default ctor.
 *
 * \param u  Scalar value.
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite(const T &u)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] = u;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for N = 2.
 *
 * \param u0 1st element.
 * \param u1 2nd element.
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite(const T &u0,
                               const T &u1)
{
    REQUIRE(N == 2);
    d_U[0] = u0;
    d_U[1] = u1;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for N = 3.
 *
 * \param u0 1st element.
 * \param u1 2nd element.
 * \param u2 3rd element.
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite(const T &u0,
                               const T &u1,
                               const T &u2)
{
    REQUIRE(N == 3);
    d_U[0] = u0;
    d_U[1] = u1;
    d_U[2] = u2;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for N = 4.
 *
 * \param u0 1st element.
 * \param u1 2nd element.
 * \param u2 3rd element.
 * \param u3 4th element.
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite(const T &u0,
                               const T &u1,
                               const T &u2,
                               const T &u3)
{
    REQUIRE(N == 4);
    d_U[0] = u0;
    d_U[1] = u1;
    d_U[2] = u2;
    d_U[3] = u3;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for N = 5.
 *
 * \param u0 1st element.
 * \param u1 2nd element.
 * \param u2 3rd element.
 * \param u3 4th element.
 * \param u4 5th element.
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite(const T &u0,
                               const T &u1,
                               const T &u2,
                               const T &u3,
                               const T &u4)
{
    REQUIRE(N == 5);
    d_U[0] = u0;
    d_U[1] = u1;
    d_U[2] = u2;
    d_U[3] = u3;
    d_U[4] = u4;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Fill in from C array.
 *
 * \param u c-style pointer to array of size N.
 *
 */
template <class T, size_t N>
void Vector_Lite<T, N>::fill(const T u[N])
{
    for (size_t i=0; i<N; ++i)
    {
        d_U[i] = u[i];
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment to another Vector_Lite.
 */
template <class T, size_t N>
Vector_Lite<T, N>& Vector_Lite<T, N>::operator=(const Vector_Lite &rhs)
{
    // return if objects are the same
    if (&rhs == this)
        return *this;

    // otherwise do member assignment
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] = rhs.d_U[i];
    }

    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment to a scalar.
 */
template <class T, size_t N>
Vector_Lite<T, N>& Vector_Lite<T, N>::operator=(const T &rhs)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] = rhs;
    }

    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Comparison to another Vector_Lite.
 */
template <class T, size_t N>
bool Vector_Lite<T, N>::operator==(const Vector_Lite<T, N> &a) const
{
    for ( size_type i = 0; i < N; ++i )
    {
        if ( d_U[i] != a.d_U[i] ) return false;
    }

    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Comparison to another Vector_Lite.
 */
template <class T, size_t N>
bool Vector_Lite<T, N>::operator<(const Vector_Lite<T, N> &a) const
{
    for ( size_type i = 0; i < N; ++i )
    {
        if ( ! ( d_U[i] <  a.d_U[i] ) ) return false;
    }

    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Element-wise not-equals.
 */
template<class T, size_t N>
inline bool operator!=(Vector_Lite<T,N> const & lhs,
                       Vector_Lite<T,N> const & rhs)
{
    return !(lhs == rhs);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Element-wise greater than.
 */
template<class T, size_t N>
inline bool operator>(Vector_Lite<T,N> const & lhs,
                      Vector_Lite<T,N> const & rhs)
{
    return rhs < lhs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Element-wise less-than-or-equal.
 */
template<class T, size_t N>
inline bool operator<=(Vector_Lite<T,N> const & lhs,
                       Vector_Lite<T,N> const & rhs)
{
    return !(lhs > rhs);
}

//---------------------------------------------------------------------------//
/*!
 * \brief  Element-wise greater-than-or-equal.
 */
template<class T, size_t N>
inline bool operator>=(Vector_Lite<T,N> const & lhs,
                       Vector_Lite<T,N> const & rhs)
{
    return !(rhs > lhs);
}

//---------------------------------------------------------------------------//
// BASIC ARITHMETIC MEMBER FUNCTIONS, VECTOR RIGHT-HAND SIDE
//---------------------------------------------------------------------------//
/*!
 * \brief Support for +=, Vector_Lite right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator+=(const Vector_Lite<T, N> &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] += a.d_U[i];
    }

    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Support for -=, Vector_Lite right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator-=(const Vector_Lite<T, N> &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] -= a.d_U[i];
    }

    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Support for *=, Vector_Lite right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator*=(const Vector_Lite<T, N> &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] *= a.d_U[i];
    }

    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Support for /=, Vector_Lite right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator/=(const Vector_Lite<T, N> &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] /= a.d_U[i];
    }

    return *this;
}

//---------------------------------------------------------------------------//
// BASIC ARITHMETIC MEMBER FUNCTIONS, SCALAR RIGHT-HAND SIDE
//---------------------------------------------------------------------------//
/*!
 * \brief Support for +=, scalar right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator+=(const T &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] += a;
    }

    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Support for -=, scalar right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator-=(const T &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] -= a;
    }

    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Support for *=, scalar right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator*=(const T &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] *= a;
    }

    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Support for /=, scalar right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator/=(const T &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] /= a;
    }

    return *this;
}

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Computes the inner product between two vectors.
 *
 * \param a 1st vector.
 * \param b 2nd vector.
 */
template <class T, size_t N>
T inner_product(const Vector_Lite<T, N> &a,
                const Vector_Lite<T, N> &b)
{
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

//---------------------------------------------------------------------------//
// GLOBAL OPERATORS
//---------------------------------------------------------------------------//
/*!
 * \brief \a a + \a b, element by element.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator+(const Vector_Lite<T, N> &a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(a) += b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a - \a b, element by element.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator-(const Vector_Lite<T, N> &a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(a) -= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a * \a b, element by element.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator*(const Vector_Lite<T, N> &a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(a) *= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a / \a b, element by element.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator/(const Vector_Lite<T, N> &a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(a) /= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a b added to all elements of \a a.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator+(const Vector_Lite<T, N> &a,
                                         const T                  b)
{
    return Vector_Lite<T, N>(a) += b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a added to all elements of \a b.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator+(const T                  a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(b) += a;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a b subracted from all elements of \a a.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator-(const Vector_Lite<T, N> &a,
                                         const T                  b)
{
    return Vector_Lite<T, N>(a) -= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a subtracted from all elements of \a b.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator-(const T                  a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(b) -= a;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a b multiplied with all elements of \a a.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator*(const Vector_Lite<T, N> &a,
                                         const T                  b)
{
    return Vector_Lite<T, N>(a) *= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a multiplied with all elements of \a b.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator*(const T                  a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(b) *= a;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a b divided into all elements of \a a.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator/(const Vector_Lite<T, N> &a,
                                         const T                  b)
{
    return Vector_Lite<T, N>(a) /= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief negates all of \a a.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator-(const Vector_Lite<T, N> &a)
{
    Vector_Lite<T, N> neg(a);

    for (size_t i = 0; i < N; ++i )
        neg(i) = -a(i);

    return neg;
}

//---------------------------------------------------------------------------//
// STREAM OPEATORS
//---------------------------------------------------------------------------//
/*!
 * \brief Write the elements of \a a to stream \a os.
 */
template <class T, size_t N>
std::ostream &operator<<(std::ostream            &os,
                         const Vector_Lite<T, N> &a)
{
    os << a(0);

    for (size_t i = 1; i < N; ++i)
    {
        os << " " << a(i);
    }

    return os;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read elements into \a a from \a is.
 */
template <class T, size_t N>
std::istream &operator>>(std::istream      &is,
                         Vector_Lite<T, N> &a)
{
    for (size_t i = 0; i < N; ++i)
    {
        is >> a(i);
    }

    return is;
}

} // end namespace profugus

#endif // utils_Vector_Lite_i_hh

//---------------------------------------------------------------------------//
//              end of utils/Vector_Lite.i.hh
//---------------------------------------------------------------------------//
