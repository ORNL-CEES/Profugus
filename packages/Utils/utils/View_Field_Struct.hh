//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Nemesis/utils/View_Field_Struct.hh
 * \author Seth R Johnson
 * \date   Mon Jan 26 20:14:12 2015
 * \brief  View_Field_Struct helper function.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Nemesis_utils_View_Field_Struct_hh
#define Nemesis_utils_View_Field_Struct_hh

#include <type_traits>
#include <cstddef>

#include "View_Field.hh"

//---------------------------------------------------------------------------//
/*!
 * \def NEMESIS_MAKE_STRUCT_VIEW
 * \brief Return a view of a member inside a view field of structs.
 *
 * As an example, suppose that we want to keep sampling data for CE
 * associated with the given energy point, but we still want to access the
 * energy grid. In that case, we can do:
 *
 * \code
    struct EnergyPoint
    {
        float  pdf;
        float  cdf;
        double energy;
    };

    std::vector<EnergyPoint> points;

    nemesis::const_View_Field<double> energies()
    {
        return NEMESIS_MAKE_STRUCT_VIEW(points, energy);
    }
 * \endcode
 *
 * Note of course that the resulting field will have stride != 1.
 *
 * \param VIEW View field of structs to extract
 * \param MEMBER Name of the member inside that struct
 */
#define NEMESIS_MAKE_STRUCT_VIEW_IMPL(VIEW, MEMBER, STRUCT_TYPE) \
    ::nemesis::detail::make_struct_view<STRUCT_TYPE, \
                                        decltype(STRUCT_TYPE::MEMBER) \
                                       >( ::nemesis::make_view(VIEW), \
                                          offsetof(STRUCT_TYPE, MEMBER))

#define NEMESIS_MAKE_STRUCT_VIEW(VIEW, MEMBER) \
    (NEMESIS_MAKE_STRUCT_VIEW_IMPL(VIEW, MEMBER, \
                                 std::decay<decltype(VIEW)>::type::value_type))


namespace nemesis
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * \brief Create a struct view from given struct and the member value offset.
 *
 * This should only be called through the NEMESIS_MAKE_STRUCT_VIEW macro.
 */
template<class S, class T>
nemesis::const_View_Field<T>
make_struct_view(nemesis::const_View_Field<S> view_struct,
                 std::size_t                  offset)
{
    static_assert(sizeof(S) % sizeof(T) == 0,
            "Given structure and data member have bad alignment!");
    Require(offset < sizeof(S));

    // Begin/end pointers for offsets
    const char* begin_ptr
        = reinterpret_cast<const char*>(view_struct.data()) + offset;
    const char* end_ptr
        = reinterpret_cast<const char*>(view_struct.data() + view_struct.size())
            + offset;

    return nemesis::const_View_Field<T>(
            reinterpret_cast<const T*>(begin_ptr),
            reinterpret_cast<const T*>(end_ptr),
            sizeof(S) / sizeof(T));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a struct view from given struct and the member value offset.
 *
 * This should only be called through the NEMESIS_MAKE_STRUCT_VIEW macro.
 */
template<class S, class T>
nemesis::View_Field<T>
make_struct_view(nemesis::View_Field<S> view_struct,
                 std::size_t            offset)
{
    static_assert(sizeof(S) % sizeof(T) == 0,
            "Given structure and data member have bad alignment!");
    Require(offset < sizeof(S));

    // Begin/end pointers for offsets
    char* begin_ptr
        = reinterpret_cast<char*>(view_struct.data()) + offset;
    char* end_ptr
        = reinterpret_cast<char*>(view_struct.data() + view_struct.size())
            + offset;

    return nemesis::View_Field<T>(
            reinterpret_cast<T*>(begin_ptr),
            reinterpret_cast<T*>(end_ptr),
            sizeof(S) / sizeof(T));

}

} // end namespace detail

//---------------------------------------------------------------------------//
} // end namespace nemesis

//---------------------------------------------------------------------------//
#endif // Nemesis_utils_View_Field_Struct_hh

//---------------------------------------------------------------------------//
// end of Nemesis/utils/View_Field_Struct.hh
//---------------------------------------------------------------------------//
