//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Box_Shape.hh
 * \author Thomas M. Evans
 * \date   Tuesday May 6 16:40:41 2014
 * \brief  Box_Shape class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Box_Shape_hh
#define mc_Box_Shape_hh

#include "Shape.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Box_Shape
 * \brief Box shape for MC geometric sampling.
 *
 * Does uniform sampling in an orthogonal box.  The sampling formulas for each
 * side of the box is:
 * \f[
   x = \Delta_{x}\xi_{1} + x_{o}
 * \f]
 * \f[
   y = \Delta_{y}\xi_{2} + y_{o}
 * \f]
 * \f[
   z = \Delta_{z}\xi_{3} + z_{o}
 * \f]
 */
//===========================================================================//

class Box_Shape : public Shape
{
  public:
    //@{
    //! Base-class typedefs.
    typedef Shape              Base;
    typedef Base::Buffer       Buffer;
    typedef Base::Space_Vector Space_Vector;
    //@}

  private:
    // >>> DATA

    // Low corner of box.
    double d_lox, d_loy, d_loz;

    // Widths of box in each direction.
    double d_Dx, d_Dy, d_Dz;

  public:
    // Constructor.
    Box_Shape(double lox, double hix, double loy,
              double hiy, double loz, double hiz);

    // Constructor from a packed buffer.
    explicit Box_Shape(Buffer &buffer);

    // >>> ACCESSORS
    //! Return the low corner of the box.
    Space_Vector low_corner() const
    {
        return Space_Vector(d_lox, d_loy, d_loz);
    }

    //! Return the high corner of the box.
    inline Space_Vector high_corner() const
    {
        return Space_Vector(d_lox + d_Dx, d_loy + d_Dy, d_loz + d_Dz);
    }

    // >>> DERIVED INTERFACE
    //! Sample a point in the shape.
    Space_Vector sample(RNG_t &rng) const
    {
        Space_Vector point(d_Dx * rng.ran() + d_lox,
                           d_Dy * rng.ran() + d_loy,
                           d_Dz * rng.ran() + d_loz);
        ENSURE(is_point_inside(point));
        return point;
    }

    // Pack the box into a buffer
    void pack(Buffer &buffer) const;

    //! Return the volume.
    double volume() const { return d_Dx * d_Dy * d_Dz; }

    // Whether a point is on or inside this shape
    bool is_point_inside(const Space_Vector& x) const;

    // Get the bounding box
    void get_bbox(
            Space_Vector& low_corner,
            Space_Vector& high_corner) const;
};

} // end namespace profugus

#endif // mc_Box_Shape_hh

//---------------------------------------------------------------------------//
//                        end of shapes/Box_Shape.hh
//---------------------------------------------------------------------------//
