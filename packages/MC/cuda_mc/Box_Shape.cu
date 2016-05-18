//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Box_Shape.cu
 * \author Steven Hamilton
 * \date   Fri Dec 13 13:30:08 2013
 * \brief  Box shape member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Box_Shape.cuh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param lox low-x boundary of box, must be \c <= the high-x boundary
 * \param hix high-x boundary of box, must be \c >= the low-x boundary
 * \param loy low-y boundary of box, must be \c <= the high-y boundary
 * \param hiy high-y boundary of box, must be \c >= the low-y boundary
 * \param loz low-z boundary of box, must be \c <= the high-z boundary
 * \param hiz high-z boundary of box, must be \c >= the low-z boundary
 */
Box_Shape::Box_Shape(double lox,
                     double hix,
                     double loy,
                     double hiy,
                     double loz,
                     double hiz)
    : d_lox(lox)
    , d_loy(loy)
    , d_loz(loz)
    , d_Dx(hix - lox)
    , d_Dy(hiy - loy)
    , d_Dz(hiz - loz)
{
    DEVICE_REQUIRE(d_Dx >= 0.0);
    DEVICE_REQUIRE(d_Dy >= 0.0);
    DEVICE_REQUIRE(d_Dz >= 0.0);
}

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                          end of Box_Shape.cu
//---------------------------------------------------------------------------//
