//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Vector_Traits.hh
 * \author Seth R Johnson
 * \date   Wed Aug 14 11:55:46 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Vector_Traits_hh
#define cuda_utils_Vector_Traits_hh

#include "Definitions.hh"

namespace profugus
{
// Declare view fields
template <typename T> class View_Field;
template <typename T> class const_View_Field;
}

namespace cuda_utils
{
// declare host/device/multi/pseudo classes
template<typename T> class Host_Vector;
template<typename Arch_Switch, typename T> class Multi_Vector;
template<typename Arch_Switch, typename T> class Device_Vector;
template<typename Arch_Switch, typename T> class Texture_Vector;

//===========================================================================//
/*!
 * \struct Vector_Traits
 * \brief  Traits class for Host/Device switching
 *
 * \tparam Arch_Switch Architecture (Host/Device)
 * \tparam Float_T  Floating point type
 * \tparam Int_T    Integer type
 */
//===========================================================================//
template<typename Arch_Switch, typename Float_T = float, class Int_T = int>
struct Vector_Traits
{
    //@{
    //! Template typedefs
    typedef Arch_Switch Arch_t;
    typedef Float_T     float_type;
    typedef Int_T       int_type;
    //@}

    //@{
    //! View field typedefs
    typedef profugus::View_Field<float_type>       View_Field_Float;
    typedef profugus::const_View_Field<float_type> const_View_Field_Float;
    typedef profugus::View_Field<int_type>         View_Field_Int;
    typedef profugus::const_View_Field<int_type>   const_View_Field_Int;
    //@}

    //@{
    //! Host vector typedefs
    typedef Host_Vector<float_type> Host_Vector_Float;
    typedef Host_Vector<int_type>   Host_Vector_Int;
    //@}

    //@{
    //! Device (or pseudodevice) vector typedefs (floats)
    typedef Device_Vector<Arch_t, float_type>  Device_Vector_Float;
    typedef Texture_Vector<Arch_t, float_type> Texture_Vector_Float;
    typedef Multi_Vector<Arch_t, float_type>   Multi_Vector_Float;
    //@}

    //@{
    //! Device (or pseudodevice) vector typedefs (ints)
    typedef Device_Vector<Arch_t, int_type>  Device_Vector_Int;
    typedef Texture_Vector<Arch_t, int_type> Texture_Vector_Int;
    typedef Multi_Vector<Arch_t, int_type>   Multi_Vector_Int;
    //@}
};

} // end namespace cuda_utils

#endif // cuda_utils_Vector_Traits_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Vector_Traits.hh
//---------------------------------------------------------------------------//
