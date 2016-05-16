//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/geometry/Bounding_Box.hh
 * \author Seth R Johnson
 * \date   Wed Jul 01 14:55:00 2015
 * \brief  Bounding_Box class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_geometry_Bounding_Box_hh
#define MC_geometry_Bounding_Box_hh

#include <iosfwd>

#include "Utils/utils/Definitions.hh"
#include "Utils/utils/Vector_Lite.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Bounding_Box
 * \brief Axis-aligned bounding box
 */
/*!
 * \example core/test/tstBounding_Box.cc
 *
 * Test of Bounding_Box.
 */
//===========================================================================//

class Bounding_Box
{
  private:
    // >>> DATA

    typedef def::Space_Vector Space_Vector;

    Space_Vector d_lower;
    Space_Vector d_upper;

  public:
    // Construct a bbox from space vectors
    Bounding_Box(Space_Vector lo, Space_Vector hi);

    // Construct a bbox.
    Bounding_Box(double lox, double hix, double loy,
                 double hiy, double loz, double hiz);

    // >>> ACCESSORS

    //! Lower bbox coordinate
    Space_Vector lower() const { return d_lower; }

    //! Upper bbox coordinate
    Space_Vector upper() const { return d_upper; }

    // Test whether a point is inside or on a bounding box.
    inline bool is_point_inside(const Space_Vector& p) const;

    // Test the intersection with another bounding box.
    inline bool intersects(const Bounding_Box& other) const;

    // Return whether we enclose bounding box B
    inline bool encloses(const Bounding_Box& other) const;

    // Calculate the center of the bounding box
    Space_Vector calc_center() const;

    // Return the union of this bounding box and another
    Bounding_Box calc_union(const Bounding_Box& other) const;
};

//---------------------------------------------------------------------------//
// FUNCTIONS
//---------------------------------------------------------------------------//

// Print to a stream
std::ostream& operator<<(std::ostream&, const Bounding_Box&);

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "Bounding_Box.i.hh"
//---------------------------------------------------------------------------//
#endif // MC_geometry_Bounding_Box_hh

//---------------------------------------------------------------------------//
// end of geometry/Bounding_Box.hh
//---------------------------------------------------------------------------//
