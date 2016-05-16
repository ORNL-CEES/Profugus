//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/geometry/Bounding_Box.i.hh
 * \author Seth R Johnson
 * \date   Wed Jul 01 14:55:00 2015
 * \brief  Bounding_Box inline method definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_geometry_Bounding_Box_i_hh
#define MC_geometry_Bounding_Box_i_hh

namespace profugus
{
//---------------------------------------------------------------------------//
/*!
 * \brief Return whether a point is inside this box.
 */
bool Bounding_Box::is_point_inside(const Space_Vector& p) const
{
    using def::X; using def::Y; using def::Z;
    return (   d_lower[X] <= p[X]
            && d_upper[X] >= p[X]
            && d_lower[Y] <= p[Y]
            && d_upper[Y] >= p[Y]
            && d_lower[Z] <= p[Z]
            && d_upper[Z] >= p[Z]
           );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Test the intersection of two bounding boxes.
 */
bool Bounding_Box::intersects(const Bounding_Box& b) const
{
    using def::X; using def::Y; using def::Z;

    if (d_upper[X] < b.d_lower[X] || d_lower[X] > b.d_upper[X]) return false;
    if (d_upper[Y] < b.d_lower[Y] || d_lower[Y] > b.d_upper[Y]) return false;
    if (d_upper[Z] < b.d_lower[Z] || d_lower[Z] > b.d_upper[Z]) return false;

    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return whether bounding box A encloses bounding box B
 */
bool Bounding_Box::encloses(const Bounding_Box& b) const
{
    using def::X; using def::Y; using def::Z;

    return (   d_lower[X] <= b.d_lower[X]
            && d_upper[X] >= b.d_upper[X]
            && d_lower[Y] <= b.d_lower[Y]
            && d_upper[Y] >= b.d_upper[Y]
            && d_lower[Z] <= b.d_lower[Z]
            && d_upper[Z] >= b.d_upper[Z]
           );
}

//---------------------------------------------------------------------------//
} // end namespace profugus

#endif // MC_geometry_Bounding_Box_i_hh

//---------------------------------------------------------------------------//
// end of geometry/Bounding_Box.i.hh
//---------------------------------------------------------------------------//
