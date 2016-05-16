//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/utils/View_Field_Fast.hh
 * \author Seth R Johnson
 * \date   Thu Mar 31 12:28:00 2016
 * \brief  View_Field_Fast class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_View_Field_Fast_hh
#define Utils_utils_View_Field_Fast_hh

namespace profugus
{

//===========================================================================//
/*!
 * \class View_Field_Fast
 * \brief Pointer-based views into data
 *
 * This allows the compiler to use extra optimizations for stride-1 data. It is
 * primiarily meant for the "range" operation.
 */
//===========================================================================//

template<typename T>
class View_Field_Fast
{
  public:
    //@{
    //! Typedefs
    typedef T             value_type;
    typedef T&            reference;
    typedef const T&      const_reference;
    typedef T*            pointer;
    typedef const T*      const_pointer;
    typedef pointer       iterator;
    typedef const_pointer const_iterator;
    typedef std::size_t   size_type;
    //@}

  private:
    // >>> DATA

    pointer d_begin;
    pointer d_end;

  public:

    // Constructor
    View_Field_Fast(pointer begin, pointer end)
        : d_begin(begin)
        , d_end(end)
    {
        /* * */
    }

    //@{
    //! Iterator access.
    iterator begin()             { return d_begin; }
    const_iterator begin() const { return d_begin; }
    iterator end()             { return d_end; }
    const_iterator end() const { return d_end; }
    //@}
};

//===========================================================================//
/*!
 * \class const_View_Field_Fast
 * \brief Pointer-based views into const data
 *
 * This allows the compiler to use extra optimizations for stride-1 data. It is
 * primiarily meant for the "range" operation.
 */
//===========================================================================//

template<typename T>
class const_View_Field_Fast
{
  public:
    //@{
    //! Typedefs
    typedef T             value_type;
    typedef T&            reference;
    typedef const T&      const_reference;
    typedef const T*      const_pointer;
    typedef const_pointer const_iterator;
    typedef std::size_t   size_type;
    //@}

  private:
    // >>> DATA

    const_pointer d_begin;
    const_pointer d_end;

  public:

    //! Constructor
    const_View_Field_Fast(const_pointer begin, const_pointer end)
        : d_begin(begin)
        , d_end(end)
    {
        /* * */
    }

    //@{
    //! Iterator access.
    const_iterator begin() const { return d_begin; }
    const_iterator end() const   { return d_end; }
    //@}
};

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
#endif // Utils_utils_View_Field_Fast_hh

//---------------------------------------------------------------------------//
// end of Utils/utils/View_Field_Fast.hh
//---------------------------------------------------------------------------//
