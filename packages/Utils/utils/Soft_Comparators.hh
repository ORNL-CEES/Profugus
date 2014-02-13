//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Soft_Comparators.hh
 * \author Gregory G. Davidson
 * \date   Wed Sep 11 21:27:43 2013
 * \brief  Defines a series of "soft" comparator functions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Soft_Comparators_hh
#define utils_Soft_Comparators_hh

#include <cmath>
#include <functional>
#include <algorithm>
#include <utility>
#include <iterator>
#include <limits>

namespace profugus
{

//===========================================================================//
// EPSILON HELPER FUNCTION
//===========================================================================//
/*!
 * \brief Return default epsilon for a given type.
 */
template<class ValType>
static ValType default_epsilon()
{
    return std::numeric_limits<ValType>::epsilon();
}

//===========================================================================//
/*!
 * \class SoftOpBase
 * \brief Base class for the Soft_Comparators.
 */
template<class ValType>
class SoftOpBase : public std::binary_function<ValType, ValType, bool>
{
  public:
    //! Useful typedef.
    typedef ValType     value_type;

  protected:
    // Stores the epsilon parameter
    value_type b_epsilon;

  public:
    //! Constructor.
    explicit SoftOpBase(value_type eps = default_epsilon<value_type>())
        : b_epsilon(eps)
    {   }

    //! Virtual destructor.
    virtual ~SoftOpBase() {  }

    //! Return the epsilon.
    value_type epsilon() const { return b_epsilon; }

    //! Pure virtual operator.
    virtual bool operator()(value_type val_1, value_type val_2) const = 0;
};

//===========================================================================//
/*!
 * \class SoftIsEqual
 * \brief Returns whether two values are equal within an optional epsilon
 *        parameter.
 */
//===========================================================================//
template<class ValType>
class SoftIsEqual : public SoftOpBase<ValType>
{
  public:
    //! Constructor.
    explicit SoftIsEqual(ValType eps = default_epsilon<ValType>())
        : SoftOpBase<ValType>(eps)
    {   }

    //! Returns whether \a val_1 and \a val_2 are equal within \a epsilon.
    bool operator()(ValType val_1, ValType val_2) const
    {
        // If the absolute difference in the values is less than epsilon, they
        // are considered equal
        if ( std::fabs(val_1 - val_2) < SoftOpBase<ValType>::b_epsilon)
        {
            return true;
        }
        // If both values are smaller than epsilon, they are considered equal
        else if(std::fabs(val_1) < SoftOpBase<ValType>::b_epsilon &&
                std::fabs(val_2) < SoftOpBase<ValType>::b_epsilon)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};

//===========================================================================//
/*!
 * \class SoftIsNotEqual
 * \brief Returns whether two values are not equal within an optional epsilon
 *        parameter.
 */
//===========================================================================//
template<class ValType>
class SoftIsNotEqual : public SoftOpBase<ValType>
{
  public:
    // Constructor
    explicit SoftIsNotEqual(ValType eps = default_epsilon<ValType>())
        : SoftOpBase<ValType>(eps)
    {   }

    //! Returns whether \a val_1 is not equal to \a val_2 within \a epsilon.
    virtual bool operator()(ValType val_1, ValType val_2) const
    {
        SoftIsEqual<ValType> soft_is_eql(SoftOpBase<ValType>::b_epsilon);
        return !soft_is_eql(val_1, val_2);
    }
};

//===========================================================================//
/*!
 * \class SoftIsLess
 * \brief Returns whether a value is less than another value including
 *        epsilon.
 */
//===========================================================================//
template<class ValType>
class SoftIsLess : public SoftOpBase<ValType>
{
  public:
    //! Constructor.
    explicit SoftIsLess(ValType eps = default_epsilon<ValType>())
        : SoftOpBase<ValType>(eps)
    {   }

    //! Return whether \a val_1 is less than \a val_2 within \a epsilon.
    virtual bool operator()(ValType val_1, ValType val_2) const
    {
        SoftIsNotEqual<ValType> is_not_eql(SoftOpBase<ValType>::b_epsilon);
        std::less<ValType> is_less;
        return is_less(val_1, val_2) && is_not_eql(val_1, val_2);
    }
};

//===========================================================================//
/*!
 * \class SoftIsGreater
 * \brief Returns whether a value is greater than another value including
 *        epsilon.
 */
//===========================================================================//
template<class ValType>
class SoftIsGreater : public SoftOpBase<ValType>
{
  public:
    //! Constructor.
    explicit SoftIsGreater(ValType eps = default_epsilon<ValType>())
        : SoftOpBase<ValType>(eps)
    {   }

    //! Return whether \a val_1 is greater than \a val_2 within \a epsilon.
    bool operator()(ValType val_1, ValType val_2) const
    {
        SoftIsNotEqual<ValType> is_not_eql(SoftOpBase<ValType>::b_epsilon);
        std::greater<ValType> is_greater;
        return is_greater(val_1, val_2) && is_not_eql(val_1, val_2);
    }
};

//===========================================================================//
/*!
 * \class SoftIsLessEqual
 * \brief Returns whether a value is less than or equal than another value
 *        including epsilon.
 */
//===========================================================================//
template<class ValType>
class SoftIsLessEqual : public SoftOpBase<ValType>
{
  public:
    //! Constructor.
    explicit SoftIsLessEqual(ValType eps = default_epsilon<ValType>())
        : SoftOpBase<ValType>(eps)
    {   }

    //! Return whether \a val_1 is less-than or equal to \a val_2 within
    //! \a epsilon.
    bool operator()(ValType val_1, ValType val_2) const
    {
        SoftIsEqual<ValType> is_eql(SoftOpBase<ValType>::b_epsilon);
        std::less<ValType> is_less;
        return is_less(val_1, val_2) || is_eql(val_1, val_2);
    }
};

//===========================================================================//
/*!
 * \class SoftIsGreaterEqual
 * \brief Returns whether a value is greater than another value including
 *        epsilon.
 */
//===========================================================================//
template<class ValType>
class SoftIsGreaterEqual : public SoftOpBase<ValType>
{
  public:
    //! Constructor.
    explicit SoftIsGreaterEqual(ValType eps = default_epsilon<ValType>())
        : SoftOpBase<ValType>(eps)
    {   }

    //! Return whether \a val_1 is greater-than or equal to \a val_2 within
    //! \a epsilon.
    bool operator()(ValType val_1, ValType val_2) const
    {
        SoftIsEqual<ValType> is_eql(SoftOpBase<ValType>::b_epsilon);
        std::greater<ValType> is_greater;
        return is_greater(val_1, val_2) || is_eql(val_1, val_2);
    }
};

//===========================================================================//
// ITERATOR TRAITS SHORTCUT CLASS
//===========================================================================//
/*!
 * \class IT_Traits
 * \brief Simple traits class for iterators.
 */
template<class InputIterator>
struct ITT
{
    typedef typename std::iterator_traits<InputIterator>::value_type vt;
    typedef std::pair<vt, vt>                                        vt_pair;
};

//===========================================================================//
// COMPARE_TRAITS CLASS
//===========================================================================//
/*!
 * \class Compare_Traits
 * \brief Simple traits class for comparators.
 */
template<class ValType>
struct Compare_Traits
{
    typedef std::less<ValType>       less_comp;
    typedef std::less_equal<ValType> less_eql_comp;
    typedef std::equal_to<ValType>   eql_comp;
};

//---------------------------------------------------------------------------//
//! Specialization for doubles.
template<>
struct Compare_Traits<double>
{
    typedef profugus::SoftIsLess<double>      less_comp;
    typedef profugus::SoftIsLessEqual<double> less_eql_comp;
    typedef profugus::SoftIsEqual<double>     eql_comp;
};

//---------------------------------------------------------------------------//
//! Specialization for floats.
template<>
struct Compare_Traits<float>
{
    typedef profugus::SoftIsLessEqual<float>  less_comp;
    typedef profugus::SoftIsLessEqual<double> less_eql_comp;
    typedef profugus::SoftIsEqual<float>      eql_comp;
};

//===========================================================================//
// HELPER FUNCTIONS
//===========================================================================//
//! Returns whether two values are equal within epsilon.
template<class ValType>
bool soft_is_equal(ValType val_1,
                   ValType val_2,
                   ValType eps = default_epsilon<ValType>())
{
    SoftIsEqual<ValType> is_equal(eps);
    return is_equal(val_1, val_2);
}

//---------------------------------------------------------------------------//
//! Returns whether two values are not equal within epsilon.
template<class ValType>
bool soft_is_not_equal(ValType val_1,
                       ValType val_2,
                       ValType eps = default_epsilon<ValType>())
{
    SoftIsNotEqual<ValType> is_not_equal(eps);
    return is_not_equal(val_1, val_2);
}

//---------------------------------------------------------------------------//
//! Returns whether val1 is less than val2 including epsilon.
template<class ValType>
bool soft_is_less(ValType val_1,
                  ValType val_2,
                  ValType eps = default_epsilon<ValType>())
{
    SoftIsLess<ValType> is_less(eps);
    return is_less(val_1, val_2);
}

//---------------------------------------------------------------------------//
//! Returns whether val1 is less-than or equal to val2 including epsilon.
template<class ValType>
bool soft_is_less_equal(ValType val_1,
                        ValType val_2,
                        ValType eps = default_epsilon<ValType>())
{
    SoftIsLessEqual<ValType> is_less_eql(eps);
    return is_less_eql(val_1, val_2);
}

//---------------------------------------------------------------------------//
//! Returns whether val1 is greater than val2 including epsilon.
template<class ValType>
inline bool soft_is_greater(ValType val_1,
                            ValType val_2,
                            ValType eps = default_epsilon<ValType>())
{
    SoftIsGreater<ValType> is_greater(eps);
    return is_greater(val_1, val_2);
}

//---------------------------------------------------------------------------//
//! Returns whether val1 is greater-than or equal to val2 including epsilon.
template<class ValType>
inline bool soft_is_greater_equal(ValType val_1,
                                  ValType val_2,
                                  ValType eps = default_epsilon<ValType>())
{
    SoftIsGreaterEqual<ValType> is_greater_eql(eps);
    return is_greater_eql(val_1, val_2);
}

//---------------------------------------------------------------------------//
//! Returns whether the values of pair1 are within the values of pair2.
template<class ValType>
inline bool soft_is_within(const std::pair<ValType, ValType> &pair1,
                           const std::pair<ValType, ValType> &pair2,
                           const ValType &eps = default_epsilon<ValType>())
{
    if( soft_is_greater_equal(pair1.first, pair2.first, eps) &&
        soft_is_less_equal(pair1.second, pair2.second, eps) )
    {
        return true;
    }
    else
    {
        return false;
    }
}

} // end namespace profugus

#endif // utils_Soft_Comparators_hh

//---------------------------------------------------------------------------//
//              end of utils/Soft_Comparators.hh
//---------------------------------------------------------------------------//
