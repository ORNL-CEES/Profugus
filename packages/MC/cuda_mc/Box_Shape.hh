//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Box_Shape.hh
 * \author Stuart Slattery
 * \brief  Box_Shape class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Box_Shape_hh
#define cuda_mc_Box_Shape_hh

#include "cuda_utils/Definitions.hh"
#include "cuda_utils/CudaMacros.hh"
#include "cuda_utils/CudaDBC.hh"

namespace cuda_profugus
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

class Box_Shape
{
  public:
    //@{
    //! Base-class typedefs.
    typedef cuda::Space_Vector Space_Vector;
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
              double hiy, double loz, double hiz)
	: d_lox(lox)
	, d_loy(loy)
	, d_loz(loz)
	, d_Dx(hix - lox)
	, d_Dy(hiy - loy)
	, d_Dz(hiz - loz)
    {
	REQUIRE(d_Dx >= 0.0);
	REQUIRE(d_Dy >= 0.0);
	REQUIRE(d_Dz >= 0.0);
    }

    // >>> ACCESSORS
    //! Return the low corner of the box.
    PROFUGUS_HOST_DEVICE_FUNCTION
    Space_Vector low_corner() const
    {
	Space_Vector r;
	r.x = d_lox;
	r.y = d_loy;
	r.z = d_loz;
	return r;
    }

    //! Return the high corner of the box.
    PROFUGUS_HOST_DEVICE_FUNCTION
    Space_Vector high_corner() const
    {
	Space_Vector r;
	r.x = d_lox + d_Dx;
	r.y = d_loy + d_Dy;
	r.z = d_loz + d_Dz;
	return r;
    }

    // >>> DERIVED INTERFACE
    //! Sample a point in the shape.
    PROFUGUS_HOST_DEVICE_FUNCTION
    Space_Vector sample(const double ran1, const double ran2, const double ran3 ) const
    {
        Space_Vector point;
	point.x = d_Dx * ran1 + d_lox;
	point.y = d_Dy * ran2 + d_loy;
	point.z = d_Dz * ran3 + d_loz;
        ENSURE(is_point_inside(point));
        return point;
    }

    //! Return the volume.
    PROFUGUS_HOST_DEVICE_FUNCTION
    double volume() const { return d_Dx * d_Dy * d_Dz; }

    // Whether a point is on or inside this shape
    PROFUGUS_HOST_DEVICE_FUNCTION
    bool is_point_inside(const Space_Vector& r) const
    {
	return (r.x >= d_lox && r.x <= d_lox + d_Dx) &&
	    (r.y >= d_loy && r.y <= d_loy + d_Dy) &&
	    (r.z >= d_loz && r.z <= d_loz + d_Dz);
    } 

    // Get the bounding box
    PROFUGUS_HOST_DEVICE_FUNCTION
    void get_bbox( Space_Vector& low_corner, Space_Vector& high_corner) const
    {
	low_corner.x = d_lox;
	low_corner.y = d_loy;
	low_corner.z = d_loz;
	high_corner.x = d_lox + d_Dx;
	high_corner.y = d_loy + d_Dy;
	high_corner.z = d_loz + d_Dz;
    }
};

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Box_Shape_hh

//---------------------------------------------------------------------------//
//                        end of cuda_mc/Box_Shape.hh
//---------------------------------------------------------------------------//
