//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/RTK_Cell.cc
 * \author Thomas M. Evans
 * \date   Tuesday April 29 15:41:12 2014
 * \brief  RTK_Cell member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <algorithm>
#include <iomanip>

#include "harness/Soft_Equivalence.hh"
#include "harness/Warnings.hh"
#include "utils/Constants.hh"
#include "utils/Vector_Functions.hh"
#include "Definitions.hh"
#include "RTK_Cell.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Empty, square pin-cell constructor.
 */
RTK_Cell::RTK_Cell(int    mod_id,
                   double pitch,
                   double height,
                   int    segments)
    : d_mod_id(mod_id)
    , d_r(0)
    , d_ids(0)
    , d_z(height)
    , d_num_shells(0)
    , d_num_regions(1)
    , d_segments(segments)
    , d_seg_faces(d_segments / 2)
    , d_num_int_faces(d_seg_faces + d_num_shells)
    , d_mod_region(d_num_regions - 1)
    , d_num_cells(d_num_regions * d_segments)
    , d_vessel(false)
{
    REQUIRE(d_z > 0.0);
    REQUIRE(d_segments == 1 || d_segments == 4);
    REQUIRE(d_mod_region == 0);

    // square cell
    d_xy[0] = pitch;
    d_xy[1] = pitch;

    d_extent[0][LO] = -d_xy[0] * 0.5;
    d_extent[0][HI] =  d_xy[0] * 0.5;
    d_extent[1][LO] = -d_xy[1] * 0.5;
    d_extent[1][HI] =  d_xy[1] * 0.5;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Empty pin-cell constructor with variable X/Y pitch.
 */
RTK_Cell::RTK_Cell(int    mod_id,
                   double dx,
                   double dy,
                   double height,
                   int    segments)
    : d_mod_id(mod_id)
    , d_r(0)
    , d_ids(0)
    , d_z(height)
    , d_num_shells(0)
    , d_num_regions(1)
    , d_segments(segments)
    , d_seg_faces(d_segments / 2)
    , d_num_int_faces(d_seg_faces + d_num_shells)
    , d_mod_region(d_num_regions - 1)
    , d_num_cells(d_num_regions * d_segments)
    , d_vessel(false)
{
    REQUIRE(d_z > 0.0);
    REQUIRE(d_segments == 1 || d_segments == 4);
    REQUIRE(d_mod_region == 0);

    // unique pitches in X/Y
    d_xy[0] = dx;
    d_xy[1] = dy;

    d_extent[0][LO] = -d_xy[0] * 0.5;
    d_extent[0][HI] =  d_xy[0] * 0.5;
    d_extent[1][LO] = -d_xy[1] * 0.5;
    d_extent[1][HI] =  d_xy[1] * 0.5;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Single shell constructor.
 */
RTK_Cell::RTK_Cell(int    fuel_id,
                   double r,
                   int    mod_id,
                   double pitch,
                   double height,
                   int    segments)
    : d_mod_id(mod_id)
    , d_r(1, r)
    , d_ids(1, fuel_id)
    , d_z(height)
    , d_num_shells(d_r.size())
    , d_num_regions(d_num_shells + 1)
    , d_segments(segments)
    , d_seg_faces(d_segments / 2)
    , d_num_int_faces(d_seg_faces + d_num_shells)
    , d_mod_region(d_num_regions - 1)
    , d_num_cells(d_num_regions * d_segments)
    , d_vessel(false)
{
    REQUIRE(d_r[0] > 0.0);
    REQUIRE(d_z > 0.0);
    REQUIRE(d_segments == 1 || d_segments == 4);
    REQUIRE(d_mod_region == 1);

    // square cell
    d_xy[0] = pitch;
    d_xy[1] = pitch;

    d_extent[0][LO] = -d_xy[0] * 0.5;
    d_extent[0][HI] =  d_xy[0] * 0.5;
    d_extent[1][LO] = -d_xy[1] * 0.5;
    d_extent[1][HI] =  d_xy[1] * 0.5;

    ENSURE(d_r[0] <= d_extent[0][HI]);
    ENSURE(d_r[0] <= d_extent[1][HI]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Multiple shell constructor.
 */
RTK_Cell::RTK_Cell(const Vec_Int &ids,
                   const Vec_Dbl &r,
                   int            mod_id,
                   double         pitch,
                   double         height,
                   int            segments)
    : d_mod_id(mod_id)
    , d_r(r)
    , d_ids(ids)
    , d_z(height)
    , d_num_shells(d_r.size())
    , d_num_regions(d_num_shells + 1)
    , d_segments(segments)
    , d_seg_faces(d_segments / 2)
    , d_num_int_faces(d_seg_faces + d_num_shells)
    , d_mod_region(d_num_regions - 1)
    , d_num_cells(d_num_regions * d_segments)
    , d_vessel(false)
{
    REQUIRE(!d_r.empty() ? d_r.front() > 0.0 : true);
    REQUIRE(d_z > 0.0);
    REQUIRE(d_segments == 1 || d_segments == 4);
    REQUIRE(d_r.size() == d_ids.size());
    REQUIRE(!d_r.empty() ? d_mod_region > 0 : d_mod_region == 0);

    // square cell
    d_xy[0] = pitch;
    d_xy[1] = pitch;

    d_extent[0][LO] = -d_xy[0] * 0.5;
    d_extent[0][HI] =  d_xy[0] * 0.5;
    d_extent[1][LO] = -d_xy[1] * 0.5;
    d_extent[1][HI] =  d_xy[1] * 0.5;

    ENSURE(!d_r.empty() ? d_r.back() <= d_extent[0][HI] : true);
    ENSURE(!d_r.empty() ? d_r.back() <= d_extent[1][HI] : true);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Multiple shell constructor with option for gap.
 *
 * \arg gap Gap_Vector organized \c (lox,hix,loy,hiy); only 1 x side and 1 y
 * side can have gap, ie. lox/hiy, lox/loy, hix/loy, hix/hiy
 */
RTK_Cell::RTK_Cell(const Vec_Int    &ids,
                   const Vec_Dbl    &r,
                   int               mod_id,
                   double            pitch,
                   double            height,
                   const Gap_Vector &gap,
                   int               segments)
    : d_mod_id(mod_id)
    , d_r(r)
    , d_ids(ids)
    , d_z(height)
    , d_num_shells(d_r.size())
    , d_num_regions(d_num_shells + 1)
    , d_segments(segments)
    , d_seg_faces(d_segments / 2)
    , d_num_int_faces(d_seg_faces + d_num_shells)
    , d_mod_region(d_num_regions - 1)
    , d_num_cells(d_num_regions * d_segments)
    , d_vessel(false)
{
    REQUIRE(d_r.front() > 0.0);
    REQUIRE(d_z > 0.0);
    REQUIRE(d_segments == 1 || d_segments == 4);
    REQUIRE(d_r.size() == d_ids.size());
    REQUIRE(d_num_shells == 1 ? d_mod_region == 1 : d_mod_region > 1);
    REQUIRE(gap[0] > 0.0 ? gap[1] == 0.0 : true);
    REQUIRE(gap[1] > 0.0 ? gap[0] == 0.0 : true);
    REQUIRE(gap[2] > 0.0 ? gap[3] == 0.0 : true);
    REQUIRE(gap[3] > 0.0 ? gap[2] == 0.0 : true);

    // square cell with gap
    d_xy[0] = pitch + gap[0] + gap[1];
    d_xy[1] = pitch + gap[2] + gap[3];

    d_extent[0][LO] = -(pitch * 0.5 + gap[0]);
    d_extent[0][HI] =   pitch * 0.5 + gap[1];
    d_extent[1][LO] = -(pitch * 0.5 + gap[2]);
    d_extent[1][HI] =   pitch * 0.5 + gap[3];

    ENSURE(d_r.back() <= pitch * 0.5);
    ENSURE(d_r.back() <= pitch * 0.5);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Empty pin-cell constructor with option for gap.
 *
 * \arg gap Gap_Vector organized \c (lox,hix,loy,hiy); only 1 x side and 1 y
 * side can have gap, ie. lox/hiy, lox/loy, hix/loy, hix/hiy
 */
RTK_Cell::RTK_Cell(int               mod_id,
                   double            pitch,
                   double            height,
                   const Gap_Vector &gap,
                   int               segments)
    : d_mod_id(mod_id)
    , d_r(0)
    , d_ids(0)
    , d_z(height)
    , d_num_shells(0)
    , d_num_regions(1)
    , d_segments(segments)
    , d_seg_faces(d_segments / 2)
    , d_num_int_faces(d_seg_faces + d_num_shells)
    , d_mod_region(d_num_regions - 1)
    , d_num_cells(d_num_regions * d_segments)
    , d_vessel(false)
{
    REQUIRE(d_z > 0.0);
    REQUIRE(d_segments == 1 || d_segments == 4);
    REQUIRE(d_mod_region == 0);
    REQUIRE(gap[0] > 0.0 ? gap[1] == 0.0 : true);
    REQUIRE(gap[1] > 0.0 ? gap[0] == 0.0 : true);
    REQUIRE(gap[2] > 0.0 ? gap[3] == 0.0 : true);
    REQUIRE(gap[3] > 0.0 ? gap[2] == 0.0 : true);

    // square cell with gap
    d_xy[0] = pitch + gap[0] + gap[1];
    d_xy[1] = pitch + gap[2] + gap[3];

    d_extent[0][LO] = -(pitch * 0.5 + gap[0]);
    d_extent[0][HI] =   pitch * 0.5 + gap[1];
    d_extent[1][LO] = -(pitch * 0.5 + gap[2]);
    d_extent[1][HI] =   pitch * 0.5 + gap[3];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Empty pin-cell constructor with variable X/Y pitch and vessel.
 *
 * An outer vessel is defined by an enclosing array.  In order to calculate
 * distance to boundary, it is necessary to know the offsets from the origin
 * of the vessel radius to the lower-left coordinates of the pin-cell, ie.
 * \verbatim
  _              _
  xo             x1 _
  |--------------/| y1
  |             / |    o is the origin of the pin cell (0,0) in local pin space
  |            /  |    O origin of the vessel cylinder
  |       o   /   |
  |          /    |
  |         /     | _
  |--------/-\----| yo
  ^           \
  |            \ R
    y_offset    \
  |              \
                  \
  |                \
  < - - - - - - - - O
       x_offset
   \endverbatim
 *
 * In order to calculate distance-to-boundary, we must transform a point in
 * the pincell (in the pincell's local coordinate system) into the vessel's
 * coordinate system.  A point in the pincell coordinate system can be defined
 * in the vessel's coordinate system using the relationships:
 * \f[
   \begin{array}{rcl}
    x &=& (\bar{x} - \bar{x}_o) + \mbox{x\_offset}\:,\\
    y &=& (\bar{y} - \bar{y}_o) + \mbox{y\_offset}\:.
   \end{array}
 * \f]
 * where \f$(\bar{x}, \bar{y})\f$ are coordinates of a point in the local
 * pin-cell system.  The origin of the pin-cell in the vessel's coordinate
 * system is thus
 * \f[
   \begin{array}{rcl}
    X_c &=& \mbox{x\_offset} - \bar{x}_o\:,\\
    Y_c &=& \mbox{y\_offset} - \bar{y}_o\:.
   \end{array}
 * \f]
 * A point in the pincell defined in the pincell's coordinate system
 * can be transformed into the vessel's coordinate system using
 * \f[
   \begin{array}{rcl}
    x &=& \bar{x} + X_c\:,\\
    y &=& \bar{y} + Y_c\:.
   \end{array}
 * \f]
 * Note, the offsets should have the correct sign, in the example above the \c
 * x_offset would be negative).
 *
 * \param R0 inner radius of the vessel
 * \param R1 outer radius of the vessel

 * \param x_offset offset between the origin of the vessel and the left-edge
 * of the pincell
 * \param y_offset offset between the origin of the vessel and the bottom-edge
 * of the pincell

 * \param vessel_id material id of the vessel
 *
 * \pre either R0, R1, or both must bisect the pincell
 */
RTK_Cell::RTK_Cell(int    mod_id,
                   double dx,
                   double dy,
                   double height,
                   double R0,
                   double R1,
                   double x_offset,
                   double y_offset,
                   int    vessel_id)
    : d_mod_id(mod_id)
    , d_r(0)
    , d_ids(0)
    , d_z(height)
    , d_num_shells(0)
    , d_num_regions(1)
    , d_segments(1)
    , d_seg_faces(d_segments / 2)
    , d_num_int_faces(d_seg_faces + d_num_shells)
    , d_mod_region(d_num_regions - 1)
    , d_num_cells(d_num_regions * d_segments)
    , d_vessel(true)
    , d_R0(-1.0)
    , d_R1(-1.0)
    , d_inner(false)
    , d_outer(false)
    , d_vessel_id(vessel_id)
{
    using def::X; using def::Y;

    REQUIRE(d_z > 0.0);
    REQUIRE(d_mod_region == 0);
    REQUIRE(R0 < R1);

    // unique pitches in X/Y
    d_xy[0] = dx;
    d_xy[1] = dy;

    d_extent[0][LO] = -d_xy[0] * 0.5;
    d_extent[0][HI] =  d_xy[0] * 0.5;
    d_extent[1][LO] = -d_xy[1] * 0.5;
    d_extent[1][HI] =  d_xy[1] * 0.5;

    // calculate offset to vessel origin from pincell origin
    d_offsets[X] = x_offset - d_extent[X][LO];
    d_offsets[Y] = y_offset - d_extent[Y][LO];

    // near and far corners relative to the origin of vessel cylinder
    double near[2], far[2];

    for (int dir = 0; dir < 2; ++dir)
    {
        near[dir] = this->l2g(d_extent[dir][LO], dir);
        far[dir]  = this->l2g(d_extent[dir][HI], dir);

        if (d_offsets[dir] < 0.0)
        {
            near[dir] = this->l2g(d_extent[dir][HI], dir);
            far[dir]  = this->l2g(d_extent[dir][LO], dir);
        }
    }

    // calculate the near and far radii bisecting the pincell
    double nearR2 = near[X]*near[X] + near[Y]*near[Y];
    double farR2  = far[X]*far[X] + far[Y]*far[Y];
    double R0_2   = R0 * R0;
    double R1_2   = R1 * R1;

    // check to see if R0 or R1 bisect the cell, R0 < R1 so if R0 > farR the
    // vessel does not bisect the cell
    VALIDATE(R0_2 < farR2,
             "R0 = " << R0 << " is greater than the far extent "
             << "of the pincell, " << std::sqrt(farR2));

    // likewise if R1 < nearR the vessel cannot bisect the cell
    VALIDATE(R1_2 > nearR2,
             "R1 = " << R1 <<  " is less than the near extent "
             << "of the pincell, " << std::sqrt(farR2));

    // now we have to check each vessel radius
    if (R0_2 > nearR2)
    {
        CHECK(R0_2 < farR2);
        d_R0    = R0;
        d_inner = true;
    }

    if (R1_2 < farR2)
    {
        CHECK(R1_2 > nearR2);
        d_R1    = R1;
        d_outer = true;
    }
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize a geometric state in the pin.
 *
 * \param r a point somewhere in this pin-cell
 * \param state geometric state that gets set
 */
void RTK_Cell::initialize(const Space_Vector &r,
                          Geo_State_t        &state) const
{
    REQUIRE(r[0] >= d_extent[0][LO]);
    REQUIRE(r[0] <= d_extent[0][HI]);
    REQUIRE(r[1] >= d_extent[1][LO]);
    REQUIRE(r[1] <= d_extent[1][HI]);
    REQUIRE(r[2] >= 0.0);
    REQUIRE(r[2] <= d_z);

    // get the current region
    state.region = region(r[0], r[1]);

    // get the current segment
    state.segment = segment(r[0], r[1]);

    // initialize to no-face (this should probably be fixed to allow particle
    // born on a face)
    state.face = Geo_State_t::NONE;

    // we are not on an exiting face
    state.exiting_face = Geo_State_t::NONE;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Track distance to boundary.
 */
void RTK_Cell::distance_to_boundary(const Space_Vector &r,
                                    const Space_Vector &omega,
                                    Geo_State_t        &state)
{
    using def::X; using def::Y; using def::Z;

    REQUIRE(!d_vessel ? state.region >= 0 && state.region < d_num_regions
             : state.region == Geo_State_t::MODERATOR ||
             state.region == Geo_State_t::VESSEL);
    REQUIRE(soft_equiv(vector_magnitude(omega), 1., 1.e-6));
    REQUIRE(omega[X]<0.0 ? r[X] >= d_extent[X][LO] : r[X] <= d_extent[X][HI]);
    REQUIRE(omega[Y]<0.0 ? r[Y] >= d_extent[Y][LO] : r[Y] <= d_extent[Y][HI]);
    REQUIRE(omega[Z]<0.0 ? r[Z] >= 0.0             : r[Z] <= d_z);

    // initialize running dist-to-boundary
    d_db                      = constants::huge;
    state.dist_to_next_region = d_db;
    state.next_segment        = state.segment;

    // >>> CHECK FOR INTERSECTIONS WITH OUTSIDE BOX

    // check radial surfaces of box
    dist_to_radial_face(X, r[X], omega[X], state);
    dist_to_radial_face(Y, r[Y], omega[Y], state);

    // check axial surface
    dist_to_axial_face(r[Z], omega[Z], state);

    // >>> CHECK FOR INTERSECTIONS WITH RADIAL SHELLS
    if (d_num_shells > 0)
    {
        calc_shell_db(r, omega, state);
    }

    // >>> CHECK FOR INTERSECTIONS WITH VESSEL
    dist_to_vessel(r, omega, state);

    // >>> CHECK FOR SEGMENT INTERSECTIONS

    // crossing a segment does not enter a different region
    if (d_segments > 1)
    {
        // initialize distance to boundary
        d_db = constants::huge;

        // check for intersection with x segment planes
        if (state.face != d_num_shells)
        {
            if (omega[X] > 0.0 && r[X] < 0.0)
            {
                d_db      = -r[X] / omega[X];
                d_face    = d_num_shells;
                d_segment = state.segment - 1;
            }
            else if (omega[X] < 0.0 && r[X] > 0.0)
            {
                d_db      = -r[X] / omega[X];
                d_face    = d_num_shells;
                d_segment = state.segment + 1;
            }

            // update distance to boundary info
            if (d_db < state.dist_to_next_region)
            {
                state.dist_to_next_region = d_db;
                state.exiting_face        = Geo_State_t::INTERNAL;
                state.next_face           = d_face;
                state.next_region         = state.region;
                state.next_segment        = d_segment;
            }
        }

        // check for intersection with y segment planes
        if (state.face != d_num_shells + 1)
        {
            if (omega[Y] > 0.0 && r[Y] < 0.0)
            {
                d_db      = -r[Y] / omega[Y];
                d_face    = d_num_shells + 1;
                d_segment = state.segment - 2;
            }
            else if (omega[Y] < 0.0 && r[Y] > 0.0)
            {
                d_db      = -r[Y] / omega[Y];
                d_face    = d_num_shells + 1;
                d_segment = state.segment + 2;
            }

            // update distance to boundary info
            if (d_db < state.dist_to_next_region)
            {
                state.dist_to_next_region = d_db;
                state.exiting_face        = Geo_State_t::INTERNAL;
                state.next_face           = d_face;
                state.next_region         = state.region;
                state.next_segment        = d_segment;
            }
        }
    }

    ENSURE(state.dist_to_next_region >= 0.0);
    ENSURE(state.exiting_face == Geo_State_t::INTERNAL ?
            state.next_region >= 0 : true);
    ENSURE(state.exiting_face == Geo_State_t::INTERNAL && !d_vessel ?
            state.next_region >= 0 && state.next_region < d_num_regions : true);
    ENSURE(state.exiting_face == Geo_State_t::INTERNAL && !d_vessel ?
            state.next_face < d_num_int_faces : true);
    ENSURE(state.next_segment >= 0 && state.next_segment < d_segments);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the geometric state during transport.
 */
void RTK_Cell::update_state(Geo_State_t &state) const
{
    // move the ray off a face
    state.face = Geo_State_t::NONE;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a particle through a surface.
 */
void RTK_Cell::cross_surface(Geo_State_t &state) const
{
    REQUIRE(state.exiting_face == Geo_State_t::INTERNAL);

    // update the current face of the particle
    state.face = state.next_face;

    // update the region of the cell that is being entered
    state.region = state.next_region;

    // update the segment
    state.segment = state.next_segment;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find the region for a point.
 */
int RTK_Cell::region(double x,
                     double y) const
{
    using def::X; using def::Y;

    REQUIRE(d_num_regions > 0);
    REQUIRE(x > d_extent[0][LO] || soft_equiv(x, d_extent[0][LO]));
    REQUIRE(x < d_extent[0][HI] || soft_equiv(x, d_extent[0][HI]));
    REQUIRE(y > d_extent[1][LO] || soft_equiv(y, d_extent[1][LO]));
    REQUIRE(y < d_extent[1][HI] || soft_equiv(y, d_extent[1][HI]));

    // if this is an empty box, return the moderator
    if (!d_num_shells)
    {
        REQUIRE(d_num_regions == 1);

        // check to see if it is in a vessel region
        if (d_vessel)
        {
            // get coordinates of point in vessel coordinate system
            double gx = this->l2g(x, X);
            double gy = this->l2g(y, Y);
            double R2 = gx*gx + gy*gy;

            if (d_inner && d_outer)
            {
                if (R2 < d_R0*d_R0)
                    return Geo_State_t::MODERATOR;
                else if (R2 > d_R1*d_R1)
                    return Geo_State_t::MODERATOR;
                else
                    return Geo_State_t::VESSEL;
            }
            else if (d_inner)
            {
                if (R2 > d_R0*d_R0)
                    return Geo_State_t::VESSEL;
                else
                    return  Geo_State_t::MODERATOR;
            }
            else if (d_outer)
            {
                if (R2 < d_R1*d_R1)
                    return Geo_State_t::VESSEL;
                else
                    return Geo_State_t::MODERATOR;
            }
            else
            {
                VALIDATE(false, "Vessel without bisection.");
            }
        }
        else
        {
            return Geo_State_t::MODERATOR;
        }
    }

    // calculate rp
    double rp2 = x*x + y*y;

    // loop through shells to find the containing shell
    for (int n = 0; n < d_num_shells; ++n)
    {
        if (rp2 <= d_r[n] * d_r[n]) return n;
    }

    // point is in the moderator region
    ENSURE(rp2 > d_r[d_num_shells-1]*d_r[d_num_shells-1]);
    return d_num_regions - 1;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get vessel data.
 *
 * \return true if vessel in cell; false otherwise
 *
 * \param R0 inner vessel radius (-1 if inner radius does not bisect vessel)
 * \param R1 outer vessel radisu (-1 if outer radius does not bisect vessel)
 * \param xc x-coordinate of cell origin relative to vessel origin
 * \param yc y-coordinate of cell origin relative to vessel origin
 */
bool RTK_Cell::vessel_data(double &R0,
                           double &R1,
                           double &xc,
                           double &yc) const
{
    // return if no vessel
    if (!d_vessel)
    {
        R0 = -1.0;
        R1 = -1.0;
        return false;
    }

    // assign the data describing the vessel
    R0 = d_R0;
    R1 = d_R1;
    xc = d_offsets[0];
    yc = d_offsets[1];
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build volumes in the pin.
 */
void RTK_Cell::build_volumes(Vec_Dbl &v,
                             int      offset) const
{
    using profugus::constants::pi;

    REQUIRE(v.size() >= offset + d_num_cells);

    // volume of each region
    Vec_Dbl vr(d_num_regions);

    // regions counter and last volume
    int     ctr = 0;
    double last = 0.0;

    // iterate through shells
    for (auto r : d_r)
    {
        // calculate the volume of this cylinder
        double cyl_v = pi * r*r * d_z;

        // assign the volume of the shell
        vr[ctr] = cyl_v - last;

        // update the last cyclinder volume
        last = cyl_v;

        // increment the ctr
        ++ctr;
    }
    CHECK(ctr == d_num_regions - 1);

    // get the volume outside the shells
    vr[ctr] = d_xy[0] * d_xy[1] * d_z - last;

    // segment-volume correction
    double seg_correction = 1.0 / static_cast<double>(d_segments);

    // assign the volumes
    for (int s = 0; s < d_segments; ++s)
    {
        for (int r = 0; r < d_num_regions; ++r)
        {
            int id = cell(r, s) + offset;
            CHECK(id < v.size());

            // correct the volume in the regions when multiple segments
            v[id] = vr[r] * seg_correction;
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Output pin cell for diagnostics.
 */
void RTK_Cell::output(std::ostream &out,
                      int           array_id,
                      int           pin_id) const
{
    using std::endl; using std::setw; using std::setprecision;
    using std::fixed; using std::left; using std::right;
    using def::I; using def::J; using def::K;

    out << "++++++++++++++++++++++++++++++++++++++++"
        << "++++++++++++++++++++++++++++++++++++++++" << endl;
    out << endl;
    out << "Pin " << setw(3) << pin_id << " from pin-cell array "
        << setw(3) << array_id << endl;
    out << endl;

    out << left;
    out << "dimensions   = " << setw(9) << setprecision(5) << fixed << d_xy[0]
        << setw(9) << setprecision(5) << fixed << d_xy[1]
        << setw(9) << setprecision(5) << fixed << d_z
        << endl;
    out << "num regions  = " << setw(4) << d_num_regions << endl;
    out << "mod-id       = " << setw(4) << d_mod_id << endl;
    if (!d_r.empty())
    {
        out << "shell radii  = ";
        for (int n = 0; n < d_r.size(); ++n)
        {
            out << setw(8) << setprecision(5) << fixed << d_r[n];
        }
        out << endl;

        out << "shell matids = ";
        for (int n = 0; n < d_r.size(); ++n)
        {
            out << setw(4) << d_ids[n];
        }
        out << endl;
    }
    out << right;

    out << endl;
    out << "++++++++++++++++++++++++++++++++++++++++"
        << "++++++++++++++++++++++++++++++++++++++++" << endl;
}

//---------------------------------------------------------------------------//
// PRIVATE INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Process distance to vessel shells.
 */
void RTK_Cell::dist_to_vessel(const Space_Vector &r,
                              const Space_Vector &omega,
                              Geo_State_t        &state)
{
    using def::X; using def::Y;

    // return if there is no vessel
    if (!d_vessel)
    {
        return;
    }

    REQUIRE(d_num_shells == 0);

    // local boolean for hitting a vessel wall
    bool hit = false;

    // check distance to first shell
    if (d_inner)
    {
        // only check if we aren't currently on the vessel face
        if (state.face != Geo_State_t::R0_VESSEL)
        {
            dist_to_shell(l2g(r[X], X), l2g(r[Y], Y), omega[X], omega[Y], d_R0,
                          Geo_State_t::R0_VESSEL);
        }

        // update the distance to boundary
        if (d_db > 0.0)
        {
            if (d_db < state.dist_to_next_region)
            {
                state.dist_to_next_region = d_db;
                state.next_face           = Geo_State_t::R0_VESSEL;
                state.exiting_face        = Geo_State_t::INTERNAL;
                hit                       = true;
            }
        }
    }

    // check distance to second shell
    if (d_outer)
    {
        // only check if we aren't currently on the vessel face
        if (state.face != Geo_State_t::R1_VESSEL)
        {
            dist_to_shell(l2g(r[X], X), l2g(r[Y], Y), omega[X], omega[Y], d_R1,
                          Geo_State_t::R1_VESSEL);
        }

        // update the distance to boundary
        if (d_db > 0.0)
        {
            if (d_db < state.dist_to_next_region)
            {
                state.dist_to_next_region = d_db;
                state.next_face           = Geo_State_t::R1_VESSEL;
                state.exiting_face        = Geo_State_t::INTERNAL;
                hit                       = true;
            }
        }
    }

    // if we are in the vessel and we hit a surface turn off the vessel;
    // otherwise we enter the vessel
    if (hit)
    {
        if (state.region == Geo_State_t::VESSEL)
        {
            state.next_region = Geo_State_t::MODERATOR;
        }
        else
        {
            state.next_region = Geo_State_t::VESSEL;
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate distance-to-internal-shell boundary.
 */
void RTK_Cell::calc_shell_db(const Space_Vector &r,
                             const Space_Vector &omega,
                             Geo_State_t        &state)
{
    REQUIRE(d_num_shells > 0);

    // if we are on a face then check both bounding faces
    if (state.face < d_num_shells)
    {
        // check to see if we hit the current shell we are on, which means
        // that we would traverse through that shells region on entrance
        if (state.region == state.face)
        {
            check_shell(r, omega, state.face, state.face, state.region + 1,
                        state.face, state);

            // if we can't hit the shell because of a glancing shot + floating
            // point error, update the region since we won't traverse the
            // shell
            if (d_db < 0.0)
            {
                state.region++;
            }
        }

        // now check for hitting shells greater than the current shell as long
        // as we aren't on the last face
        if (state.face < d_num_shells - 1)
        {
            check_shell(r, omega, state.face + 1, Geo_State_t::NONE,
                        state.region + 1, state.face + 1, state);
        }

        // now check for hitting shells less than the current shell as long as
        // we aren't on the first face
        if (state.face > 0)
        {
            check_shell(r, omega, state.face - 1, Geo_State_t::NONE,
                        state.region - 1, state.face - 1, state);
        }
    }

    // otherwise, find out where we are and check the bounding shells
    else
    {
        // check for hitting lowest shell
        if (state.region == 0)
        {
            // we can only hit the lowest shell
            check_shell(r, omega, 0, Geo_State_t::NONE, 1, 0, state);
        }

        // check for hitting highest shell
        else if (state.region == d_mod_region)
        {
            // we can only hit the outer shell
            check_shell(r, omega, d_num_shells - 1, Geo_State_t::NONE,
                        d_num_shells - 1, d_num_shells - 1, state);
        }

        // otherwise we are between shells
        else
        {
            CHECK(state.region - 1 >= 0);

            // check hitting lower shell
            check_shell(r, omega, state.region - 1, Geo_State_t::NONE,
                        state.region - 1, state.region - 1, state);

            // check hitting higher shell
            check_shell(r, omega, state.region, Geo_State_t::NONE,
                        state.region + 1, state.region, state);
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check a shell for distance to boundary.
 */
void RTK_Cell::check_shell(const Space_Vector &r,
                           const Space_Vector &omega,
                           int                 shell,
                           int                 face,
                           int                 next_region,
                           int                 next_face,
                           Geo_State_t        &state)
{
    using def::X; using def::Y; using def::Z;

    REQUIRE(shell >= 0 && shell < d_num_shells);

    // calculate the distance to the requested shell
    dist_to_shell(r[X], r[Y], omega[X], omega[Y], d_r[shell], face);

    // check the distance to boundary
    //    a) if it intersects the shell, and
    //    b) if it is the smallest distance
    if (d_db > 0.0)
    {
        if (d_db < state.dist_to_next_region)
        {
            state.dist_to_next_region = d_db;
            state.next_region         = next_region;
            state.next_face           = next_face;
            state.exiting_face        = Geo_State_t::INTERNAL;
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Distance to a shell.
 *
 * d_db is set negative if there is no intersection with the shell.
 */
void RTK_Cell::dist_to_shell(double x,
                             double y,
                             double omega_x,
                             double omega_y,
                             double r,
                             int    face)
{
    // initialize distance to boundary
    d_db = -1.0;

    // calculate terms in the quadratic
    double a = omega_x * omega_x + omega_y * omega_y;
    double b = 2.0 * (x * omega_x + y * omega_y);
    double c = x * x + y * y - r * r;

    // discriminant of quadratic
    double discriminant = b * b - 4.0 * a * c;

    // check for intersection with the surface anywhere along the line which
    // will be true if the discriminant > 0.0 (of course that doesn't mean the
    // ray will intersect the surface, just the line that the ray is on
    // intersects the surface)
    if (discriminant >= 0.0)
    {
        // calculate the sqrt of the discriminant and denominator once
        double sqr_root    = std::sqrt(discriminant);
        double denominator = 0.5 / a;

        // calculate the roots of the equation
        double d1 = (-b + sqr_root) * denominator;
        double d2 = (-b - sqr_root) * denominator;

        // determine d, if both d1 and d2 < 0 then the ray does not intersect
        // the surface
        if (d1 < 0.0)
            d_db = d2;
        else if (d2 < 0.0)
            d_db = d1;
        else if (face < d_num_shells)
            d_db = std::max(d1, d2);
        else
            d_db = std::min(d1, d2);
    }
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of RTK_Cell.cc
//---------------------------------------------------------------------------//
