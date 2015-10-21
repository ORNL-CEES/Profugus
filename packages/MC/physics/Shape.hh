//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Shape.hh
 * \author Thomas M. Evans
 * \date   Tuesday May 6 16:39:24 2014
 * \brief  Shape class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Shape_hh
#define mc_Shape_hh

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "utils/Constants.hh"
#include "utils/Vector_Lite.hh"
#include "rng/RNG.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Shape
 * \brief Defines a base class for simple geometric primitives that can be
 * sampled in MC applications.
 *
 * The Shape class defines an interface for geometric primitives (cubes,
 * planes, spheres, etc.) that can be sampled in MC applications.  The basic
 * use-case is:
 * \code
   RNG    rng;
   Box    box; // provides uniform sampling in the box
   Shape &shape = &box;

   // sample a point (uniformly) in the box
   Space_Vector point = box.sample(rng);
 * \endcode
 * The derived shape-type determines the sampling.  For example, a Tilted_Box
   derived class could be defined that samples a box using an internally
   defined shape function.  Each derived shape controls its own sampling
   strategy.
 */
//===========================================================================//

class Shape
{
  public:
    //! Space-vector typedef.
    typedef def::Space_Vector Space_Vector;

    //! Random number generator.
    typedef RNG RNG_t;

    //! Buffer typedef.
    typedef std::vector<char> Buffer;

  public:
    // Constructor.
    Shape() {  }

    // Virtual destructor.
    virtual ~Shape() {  }

    // >>> PUBLIC INTERFACE
    //! Pack the shape into a buffer
    virtual void pack(Buffer& buffer) const = 0;

    //! Sample a point in the shape.
    virtual Space_Vector sample(
	const double r1, const double r2, const double r3 ) const = 0;

    //! Volume of shape.
    virtual double volume() const = 0;

    //! Whether a point is on the surface of or inside this shape
    virtual bool is_point_inside(const Space_Vector& x) const = 0;

    //! Get the bounding box (closed set)
    virtual void get_bbox(
            Space_Vector& low_corner,
            Space_Vector& high_corner) const = 0;
};

} // end namespace profugus

#endif // mc_Shape_hh

//---------------------------------------------------------------------------//
//                        end of mc/Shape.hh
//---------------------------------------------------------------------------//
