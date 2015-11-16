//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Nemesis/utils/View_Field_Vector_Lite.hh
 * \author Seth R Johnson
 * \date   Mon Mar 16 15:35:15 2015
 * \brief  View_Field adapters for Vector_Lite
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Nemesis_utils_View_Field_Vector_Lite_hh
#define Nemesis_utils_View_Field_Vector_Lite_hh

#include "View_Field.hh"
#include "Vector_Lite.hh"

namespace nemesis
{

//---------------------------------------------------------------------------//
/*!
 * \brief Wraps a Vector_Lite in a const_View_Field.
 */
template<typename T, size_t N>
inline View_Field<T> make_view(Vector_Lite<T,N> &vec)
{
    // Empty vector gets an empty view field (data() may be undefined)
    if (vec.empty())
        return View_Field<T>();

    return View_Field<T>(vec.data(), vec.data() + vec.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wraps a Vector_Lite in a const_View_Field.
 */
template<typename T, size_t N>
inline const_View_Field<T> make_view(const Vector_Lite<T,N> &vec)
{
    // Empty vector gets an empty view field (data() may be undefined)
    if (vec.empty())
        return const_View_Field<T>();

    return const_View_Field<T>(vec.data(), vec.data() + vec.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wraps a Vector_Lite in a const_View_Field.
 */
template<typename T, size_t N>
inline const_View_Field<T> make_const_view(const Vector_Lite<T,N> &vec)
{
    return make_view<T,N>(vec);
}

//---------------------------------------------------------------------------//
} // end namespace nemesis

#endif // Nemesis_utils_View_Field_Vector_Lite_hh

//---------------------------------------------------------------------------//
// end of Nemesis/utils/View_Field_Vector_Lite.hh
//---------------------------------------------------------------------------//
