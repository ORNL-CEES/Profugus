//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/utils/View_Field_Teuchos.hh
 * \author Seth R Johnson
 * \date   Sat Apr 16 14:29:03 2016
 * \brief  View_Field_Teuchos class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_View_Field_Teuchos_hh
#define Utils_utils_View_Field_Teuchos_hh

#include <Teuchos_Array.hpp>

namespace profugus
{
//---------------------------------------------------------------------------//
/*!
 * \brief Wraps a Teuchos::Array in a const_View_Field.
 */
template<typename T>
View_Field<T> make_view(Teuchos::Array<T> &vec)
{
    // Empty vector gets an empty view field (data() may be undefined)
    if (vec.empty())
        return View_Field<T>();

    return View_Field<T>(&vec.front(), &vec.back() + 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wraps a Teuchos::Array in a const_View_Field.
 */
template<typename T>
const_View_Field<T> make_view(const Teuchos::Array<T> &vec)
{
    // Empty vector gets an empty view field (data() may be undefined)
    if (vec.empty())
        return const_View_Field<T>();

    return const_View_Field<T>(&vec.front(), &vec.back() + 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wraps a Teuchos::Array in a const_View_Field.
 */
template<typename T>
const_View_Field<T> make_const_view(const Teuchos::Array<T> &vec)
{
    return make_view<T>(vec);
}

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
#endif // Utils_utils_View_Field_Teuchos_hh

//---------------------------------------------------------------------------//
// end of Utils/utils/View_Field_Teuchos.hh
//---------------------------------------------------------------------------//
