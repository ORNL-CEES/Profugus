//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/Vector_Lite.i.hh
 * \author Thomas M. Evans
 * \date   Thu Jan  3 11:39:29 2008
 * \brief  Member definitions of class Vector_Lite.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Vector_Lite_i_hh
#define Utils_utils_Vector_Lite_i_hh

namespace profugus
{

//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor
 *
 * Initializes all values to T().
 *
 * This separate constructor is necessary for building types that have no (or
 * a 'deleted') copy constructor.
 *
 * \param u  Scalar value.
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite()
{
    for (size_type i = 0; i < N; ++i)
    {
        d_U[i] = T();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor based on a scalar value.
 *
 * Initializes all values to \a u. *
 * \param u  Scalar value.
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite(const T &u)
{
    for (size_type i = 0; i < N; ++i)
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
    static_assert(N == 2, "Invalid constructor");
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
    static_assert(N == 3, "Invalid constructor");
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
    static_assert(N == 4, "Invalid constructor");
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
    static_assert(N == 5, "Invalid constructor");
    d_U[0] = u0;
    d_U[1] = u1;
    d_U[2] = u2;
    d_U[3] = u3;
    d_U[4] = u4;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initializer list construction.
 *
 * This constructor allows initialization using C++-11 initializer lists:
 * \code
   Vector_Lite<int, 6> x = {1,2,3,4,5,6};
   \endcode
 * The following construction will work as well:
 * \code
   Vector_Lite<int, 3> x = {1,2};
   x[0] == 1;  // true
   x[1] == 2;  // true
   x[2] == 0;  // true
 * \endcode
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite(std::initializer_list<T> list)
{
    REQUIRE(list.size() <= N);
    std::fill(d_U, d_U + N, T());
    std::copy(list.begin(), list.end(), d_U);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment to a scalar.
 */
template <class T, size_t N>
Vector_Lite<T, N>& Vector_Lite<T, N>::operator=(const T &rhs)
{
    for (size_type i = 0; i < N; ++i)
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
    for (size_type i = 0; i < N; ++i)
    {
        if (d_U[i] != a.d_U[i])
            return false;
    }

    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Lexicographic comparison.
 */
template <class T, size_t N>
bool Vector_Lite<T, N>::operator<(const Vector_Lite<T, N> &a) const
{
    for (size_type i = 0; i != N; ++i)
    {
        if (d_U[i] < a.d_U[i])
            return true;
        if (d_U[i] > a.d_U[i])
            return false;
    }

    // exactly equal
    return false;
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
 * \brief Lexicographic greater than.
 */
template<class T, size_t N>
inline bool operator>(Vector_Lite<T,N> const & lhs,
                      Vector_Lite<T,N> const & rhs)
{
    return rhs < lhs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Lexicographic less-than-or-equal.
 */
template<class T, size_t N>
inline bool operator<=(Vector_Lite<T,N> const & lhs,
                       Vector_Lite<T,N> const & rhs)
{
    return !(lhs > rhs);
}

//---------------------------------------------------------------------------//
/*!
 * \brief  Lexicographic greater-than-or-equal.
 */
template<class T, size_t N>
inline bool operator>=(Vector_Lite<T,N> const & lhs,
                       Vector_Lite<T,N> const & rhs)
{
    return !(rhs > lhs);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Element-wise less-than.
 *
 * Use this function for bounds checking.
 */
template <class T, size_t N>
bool Vector_Lite<T, N>::all_lt(const Vector_Lite<T, N> &a) const
{
    for (size_type i = 0; i != N; ++i)
    {
        if (!(d_U[i] < a.d_U[i]))
            return false;
    }
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Element-wise greater-than.
 *
 * Use this function for bounds checking.
 */
template <class T, size_t N>
bool Vector_Lite<T, N>::all_gt(const Vector_Lite<T, N> &a) const
{
    return a.all_lt(*this);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Element-wise less-than-or-equal.
 *
 * Use this function for bounds checking.
 */
template <class T, size_t N>
bool Vector_Lite<T, N>::all_le(const Vector_Lite<T, N> &a) const
{
    for (size_type i = 0; i != N; ++i)
    {
        if (!(d_U[i] <= a.d_U[i]))
            return false;
    }
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Element-wise greater-than-or-equal.
 *
 * Use this function for bounds checking.
 */
template <class T, size_t N>
bool Vector_Lite<T, N>::all_ge(const Vector_Lite<T, N> &a) const
{
    return a.all_le(*this);
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
    for (size_type i = 0; i < N; ++i)
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
    for (size_type i = 0; i < N; ++i)
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
    for (size_type i = 0; i < N; ++i)
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
    for (size_type i = 0; i < N; ++i)
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
template <class T2>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator+=(const T2 &a)
{
    for (size_type i = 0; i < N; ++i)
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
template <class T2>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator-=(const T2 &a)
{
    for (size_type i = 0; i < N; ++i)
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
template <class T2>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator*=(const T2 &a)
{
    for (size_type i = 0; i < N; ++i)
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
template <class T2>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator/=(const T2 &a)
{
    for (size_type i = 0; i < N; ++i)
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
template <class T, class T2, size_t N>
inline const Vector_Lite<T, N> operator+(const Vector_Lite<T, N> &a,
                                         const T2                 b)
{
    return Vector_Lite<T, N>(a) += b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a added to all elements of \a b.
 */
template <class T, class T2, size_t N>
inline const Vector_Lite<T, N> operator+(const T2                 a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(b) += a;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a b subracted from all elements of \a a.
 */
template <class T, class T2, size_t N>
inline const Vector_Lite<T, N> operator-(const Vector_Lite<T, N> &a,
                                         const T2                 b)
{
    return Vector_Lite<T, N>(a) -= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a subtracted from all elements of \a b.
 */
template <class T, class T2, size_t N>
inline const Vector_Lite<T, N> operator-(const T2                 a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(b) -= a;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a b multiplied with all elements of \a a.
 */
template <class T, class T2, size_t N>
inline const Vector_Lite<T, N> operator*(const Vector_Lite<T, N> &a,
                                         const T2                 b)
{
    return Vector_Lite<T, N>(a) *= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a multiplied with all elements of \a b.
 */
template <class T, class T2, size_t N>
inline const Vector_Lite<T, N> operator*(const T2                 a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(b) *= a;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a b divided into all elements of \a a.
 */
template <class T, class T2, size_t N>
inline const Vector_Lite<T, N> operator/(const Vector_Lite<T, N> &a,
                                         const T2                 b)
{
    REQUIRE(b != T2(0));

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

    for (size_t i = 0; i < N; ++i)
        neg[i] = -a[i];

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
    os << a[0];

    for (size_t i = 1; i < N; ++i)
    {
        os << " " << a[i];
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
        is >> a[i];
    }

    return is;
}

} // end namespace profugus

#endif // Utils_utils_Vector_Lite_i_hh

//---------------------------------------------------------------------------//
//              end of utils/Vector_Lite.i.hh
//---------------------------------------------------------------------------//
