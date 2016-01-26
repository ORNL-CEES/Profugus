//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   geometry/Bounding_Box.cc
 * \author Seth R Johnson
 * \date   Wed Jul 01 14:55:00 2015
 * \brief  Bounding_Box class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Bounding_Box.hh"

#include <ostream>
#include <algorithm>

using def::X; using def::Y; using def::Z;

namespace profugus
{
//---------------------------------------------------------------------------//
/*!
 * \brief Construct a bbox.
 */
Bounding_Box::Bounding_Box(Space_Vector lo, Space_Vector hi)
    : d_lower(lo)
    , d_upper(hi)
{
    REQUIRE(lo[X] <= hi[X]);
    REQUIRE(lo[Y] <= hi[Y]);
    REQUIRE(lo[Z] <= hi[Z]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct a bbox from points.
 */
Bounding_Box::Bounding_Box(double lox, double hix, double loy,
                           double hiy, double loz, double hiz)
    : d_lower(lox, loy, loz)
    , d_upper(hix, hiy, hiz)
{
    REQUIRE(lox <= hix);
    REQUIRE(loy <= hiy);
    REQUIRE(loz <= hiz);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the center point.
 */
auto Bounding_Box::calc_center() const -> Space_Vector
{
    Space_Vector result((d_lower[X] + d_upper[X]) / 2.,
                        (d_lower[Y] + d_upper[Y]) / 2.,
                        (d_lower[Z] + d_upper[Z]) / 2.);
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the union of this and another bounding box
 */
Bounding_Box Bounding_Box::calc_union(const Bounding_Box& other) const
{
    Bounding_Box result(*this);
    for (int ax = 0; ax < 3; ++ax)
    {
        result.d_lower[ax] = std::min(d_lower[ax], other.d_lower[ax]);
        result.d_upper[ax] = std::max(d_upper[ax], other.d_upper[ax]);
    }
    return result;
}

/*----------------------------------------------------------------------------*/
// FREE FUNCTIONS
/*----------------------------------------------------------------------------*/
/*!
 * \brief Output to screen.
 */
std::ostream& operator<<(std::ostream& os, const Bounding_Box& bb)
{
    os << "{" << bb.lower() << " to " << bb.upper() << "}";
    return os;
}

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
// end of geometry/Bounding_Box.cc
//---------------------------------------------------------------------------//
