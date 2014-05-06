//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Box_Shape.cc
 * \author Thomas M. Evans
 * \date   Fri Dec 13 13:30:08 2013
 * \brief  Box shape member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "utils/Packing_Utils.hh"
#include "Box_Shape.hh"

namespace profugus
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
    Require (d_Dx >= 0.0);
    Require (d_Dy >= 0.0);
    Require (d_Dz >= 0.0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.  Builds from a packed buffer.
 */
Box_Shape::Box_Shape(Buffer &buffer)
{
    // Build the unpacker
    profugus::Unpacker u;
    u.set_buffer(buffer.size(), &buffer[0]);

    // Get the data
    u >> d_lox >> d_loy >> d_loz >> d_Dx >> d_Dy >> d_Dz;

    Require (d_Dx >= 0.0);
    Require (d_Dy >= 0.0);
    Require (d_Dz >= 0.0);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Pack the the box into \a buffer.
 */
void Box_Shape::pack(Buffer &buffer) const
{
    // Build the packer
    profugus::Packer p;
    p.compute_buffer_size_mode();

    // Compute the size of the stream
    p << d_lox << d_loy << d_loz << d_Dx << d_Dy << d_Dz;

    // Resize the buffer
    buffer.resize(p.size());

    // Set the buffer
    p.set_buffer(buffer.size(), &buffer[0]);

    // Pack the data
    p << d_lox << d_loy << d_loz << d_Dx << d_Dy << d_Dz;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns whether a point is on or inside this box.
 */
bool Box_Shape::is_point_inside(const Space_Vector& x) const
{
    return (x[0] >= d_lox && x[0] <= d_lox + d_Dx) &&
           (x[1] >= d_loy && x[1] <= d_loy + d_Dy) &&
           (x[2] >= d_loz && x[2] <= d_loz + d_Dz);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the bounding box.
 */
void Box_Shape::get_bbox(Space_Vector& low_corner,
                         Space_Vector& high_corner) const
{
    low_corner[0] = d_lox;
    low_corner[1] = d_loy;
    low_corner[2] = d_loz;
    high_corner[0] = d_lox + d_Dx;
    high_corner[1] = d_loy + d_Dy;
    high_corner[2] = d_loz + d_Dz;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                          end of Box_Shape.cc
//---------------------------------------------------------------------------//
