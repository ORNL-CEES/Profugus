//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_Cell.i.cuh
 * \author Tom Evans
 * \date   Mon Nov 28 12:33:05 2016
 * \brief  RTK_Cell kernel declarations.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_RTK_Cell_i_cuh
#define MC_cuda_rtk_RTK_Cell_i_cuh

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Find the cell for a point.
 *
 * \param region region id
 * \param segment segment id
 */
__device__
int RTK_Cell::cell(
    int region,
    int segment) const
{
    DEVICE_REQUIRE(d_num_regions > 0);
    DEVICE_REQUIRE(region < d_num_regions);
    DEVICE_REQUIRE(segment < d_segments);
    DEVICE_REQUIRE(d_segments == 1 || d_segments == 4);

    DEVICE_ENSURE(region + d_num_regions * segment < d_num_cells);
    return region + d_num_regions * segment;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find the cell for a point.
 */
__device__
inline int RTK_Cell::cellid(const Geo_State_Vector_t& state_vector,
                            int                       pid) const
{
    return cell(state_vector.region(pid),state_vector.segment(pid));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the matid for a region.
 */
__device__
int RTK_Cell::matid(
    int region) const
{
    DEVICE_REQUIRE(d_num_regions > 0);
    DEVICE_REQUIRE(!d_vessel ?
                   region >= 0 && region < d_num_regions :
                   region == Geo_State_t::MODERATOR ||
                   region == Geo_State_t::VESSEL);

    // return if region is in the shells
    if (region < d_num_shells)
        return d_ids[region];

    // return the vessel-id if the region is vessel
    if (region == Geo_State_t::VESSEL)
    {
        DEVICE_CHECK(d_vessel);
        return d_vessel_id;
    }

    // we are in the moderator
    return d_mod_id;
}

//---------------------------------------------------------------------------//
/*
 * \brief Get the extents in the current reference frame
 */
__device__
void RTK_Cell::get_extents(
    Space_Vector &lower,
    Space_Vector &upper) const
{
    using def::X; using def::Y; using def::Z;
    lower[X] = d_extent[X][LO]; upper[X] = d_extent[X][HI];
    lower[Y] = d_extent[Y][LO]; upper[Y] = d_extent[Y][HI];
    lower[Z] = 0.0            ; upper[Z] = d_z            ;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize a geometric state in the pin.
 *
 * \param r a point somewhere in this pin-cell
 * \param state geometric state that gets set
 */
__device__
void RTK_Cell::initialize(
    const Space_Vector &r,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    DEVICE_REQUIRE(r[0] >= d_extent[0][LO]);
    DEVICE_REQUIRE(r[0] <= d_extent[0][HI]);
    DEVICE_REQUIRE(r[1] >= d_extent[1][LO]);
    DEVICE_REQUIRE(r[1] <= d_extent[1][HI]);
    DEVICE_REQUIRE(r[2] >= 0.0);
    DEVICE_REQUIRE(r[2] <= d_z);

    // get the current region
    state_vector.region(pid) = region(r[0], r[1]);

    // get the current segment
    state_vector.segment(pid) = segment(r[0], r[1]);

    // initialize to no-face (this should probably be fixed to allow particle
    // born on a face)
    state_vector.face(pid) = Geo_State_t::NONE;

    // we are not on an exiting face
    state_vector.exiting_face(pid) = Geo_State_t::NONE;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Track distance to boundary.
 */
__device__
void RTK_Cell::distance_to_boundary(
    const Space_Vector &r,
    const Space_Vector &omega,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    using def::X; using def::Y; using def::Z;
    using cuda::utility::vector_magnitude;
    using cuda::utility::soft_equiv;

    DEVICE_REQUIRE(!d_vessel ? state_vector.region(pid) >= 0 &&
                               state_vector.region(pid) < d_num_regions
                   : state_vector.region(pid) == Geo_State_t::MODERATOR ||
                     state_vector.region(pid) == Geo_State_t::VESSEL);
    DEVICE_REQUIRE(soft_equiv(vector_magnitude(omega), 1., 1.e-6));
    DEVICE_REQUIRE(omega[X]<0.0            ?
                   r[X] >= d_extent[X][LO] :
                   r[X] <= d_extent[X][HI]);
    DEVICE_REQUIRE(omega[Y]<0.0            ?
                   r[Y] >= d_extent[Y][LO] :
                   r[Y] <= d_extent[Y][HI]);
    DEVICE_REQUIRE(omega[Z]<0.0 ?
                   r[Z] >= 0.0  :
                   r[Z] <= d_z);

    // initialize running dist-to-boundary
    double db                 = profugus::constants::huge;
    int    face               = 0;
    int    segment            = 0;
    state_vector.dist_to_next_region(pid) = db;
    state_vector.next_segment(pid)        = state_vector.segment(pid);

    // >>> CHECK FOR INTERSECTIONS WITH OUTSIDE BOX

    // check radial surfaces of box
    dist_to_radial_face(X, r[X], omega[X], state_vector, pid);
    dist_to_radial_face(Y, r[Y], omega[Y], state_vector, pid);

    // check axial surface
    dist_to_axial_face(r[Z], omega[Z], state_vector, pid);

    // >>> CHECK FOR INTERSECTIONS WITH RADIAL SHELLS
    if (d_num_shells > 0)
    {
        calc_shell_db(r, omega, state_vector, pid);
    }

    // >>> CHECK FOR INTERSECTIONS WITH VESSEL
    dist_to_vessel(r, omega, state_vector, pid);

    // >>> CHECK FOR SEGMENT INTERSECTIONS

    // crossing a segment does not enter a different region
    if (d_segments > 1)
    {
        // initialize distance to boundary
        db = profugus::constants::huge;

        // check for intersection with x segment planes
        if (state_vector.face(pid) != d_num_shells)
        {
            if (omega[X] > 0.0 && r[X] < 0.0)
            {
                db      = -r[X] / omega[X];
                face    = d_num_shells;
                segment = state_vector.segment(pid) - 1;
            }
            else if (omega[X] < 0.0 && r[X] > 0.0)
            {
                db      = -r[X] / omega[X];
                face    = d_num_shells;
                segment = state_vector.segment(pid) + 1;
            }

            // update distance to boundary info
            if (db < state_vector.dist_to_next_region(pid))
            {
                state_vector.dist_to_next_region(pid) = db;
                state_vector.exiting_face(pid)        = Geo_State_t::INTERNAL;
                state_vector.next_face(pid)           = face;
                state_vector.next_region(pid)         =
                    state_vector.region(pid);
                state_vector.next_segment(pid)        = segment;
            }
        }

        // check for intersection with y segment planes
        if (state_vector.face(pid) != d_num_shells + 1)
        {
            if (omega[Y] > 0.0 && r[Y] < 0.0)
            {
                db      = -r[Y] / omega[Y];
                face    = d_num_shells + 1;
                segment = state_vector.segment(pid) - 2;
            }
            else if (omega[Y] < 0.0 && r[Y] > 0.0)
            {
                db      = -r[Y] / omega[Y];
                face    = d_num_shells + 1;
                segment = state_vector.segment(pid) + 2;
            }

            // update distance to boundary info
            if (db < state_vector.dist_to_next_region(pid))
            {
                state_vector.dist_to_next_region(pid) = db;
                state_vector.exiting_face(pid)        = Geo_State_t::INTERNAL;
                state_vector.next_face(pid)           = face;
                state_vector.next_region(pid)         =
                    state_vector.region(pid);
                state_vector.next_segment(pid)        = segment;
            }
        }
    }

    DEVICE_ENSURE(state_vector.dist_to_next_region(pid) >= 0.0);
    DEVICE_ENSURE(state_vector.exiting_face(pid) == Geo_State_t::INTERNAL ?
                  state_vector.next_region(pid) >= 0 : true);
    DEVICE_ENSURE(state_vector.exiting_face(pid) == Geo_State_t::INTERNAL && 
                  !d_vessel ?
                  state_vector.next_region(pid) >= 0 &&
                  state_vector.next_region(pid) < d_num_regions :
                  true);
    DEVICE_ENSURE(state_vector.exiting_face(pid) == Geo_State_t::INTERNAL &&
                  !d_vessel ?
                  state_vector.next_face(pid) < d_num_int_faces : true);
    DEVICE_ENSURE(state_vector.next_segment(pid) >= 0 &&
                  state_vector.next_segment(pid) < d_segments);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the geometric state during transport.
 */
__device__
void RTK_Cell::update_state(
    Geo_State_Vector_t &state_vector, int pid) const
{
    // move the ray off a face
    state_vector.face(pid) = Geo_State_t::NONE;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a particle through a surface.
 */
__device__
void RTK_Cell::cross_surface(
    Geo_State_Vector_t &state_vector, int pid) const
{
    DEVICE_REQUIRE(state_vector.exiting_face(pid) == Geo_State_t::INTERNAL);

    // update the current face of the particle
    state_vector.face(pid) = state_vector.next_face(pid);

    // update the region of the cell that is being entered
    state_vector.region(pid) = state_vector.next_region(pid);

    // update the segment
    state_vector.segment(pid) = state_vector.next_segment(pid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find the region for a point.
 */
__device__
int RTK_Cell::region(
    double x,
    double y) const
{
    using def::X; using def::Y;
    using cuda::utility::soft_equiv;

    DEVICE_REQUIRE(d_num_regions > 0);
    DEVICE_REQUIRE(x > d_extent[0][LO] || soft_equiv(x, d_extent[0][LO]));
    DEVICE_REQUIRE(x < d_extent[0][HI] || soft_equiv(x, d_extent[0][HI]));
    DEVICE_REQUIRE(y > d_extent[1][LO] || soft_equiv(y, d_extent[1][LO]));
    DEVICE_REQUIRE(y < d_extent[1][HI] || soft_equiv(y, d_extent[1][HI]));

    // if this is an empty box, return the moderator
    if (!d_num_shells)
    {
        DEVICE_REQUIRE(d_num_regions == 1);

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
                DEVICE_INSIST(false, "Vessel without bisection.");
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
    DEVICE_ENSURE(rp2 > d_r[d_num_shells-1]*d_r[d_num_shells-1]);
    return d_num_regions - 1;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find the segment for a point
 */
__device__
int RTK_Cell::segment(
    double x,
    double y) const
{
    using def::X; using def::Y;

    DEVICE_REQUIRE(d_segments == 1 || d_segments == 4);

    // search segments
    if (d_segments == 4)
    {
        if (y < 0.0)
        {
            if (x < 0.0)
                return 3;
            else
                return 2;
        }
        else
        {
            if (x < 0.0)
                return 1;
        }
    }

    // we are in segment 0 if we get here
    return 0;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate distance-to-internal-shell boundary.
 */
__device__
void RTK_Cell::calc_shell_db(
    const Space_Vector &r,
    const Space_Vector &omega,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    DEVICE_REQUIRE(d_num_shells > 0);

    double db = profugus::constants::huge;

    // if we are on a face then check both bounding faces
    if (state_vector.face(pid) < d_num_shells)
    {
        // check to see if we hit the current shell we are on, which means
        // that we would traverse through that shells region on entrance
        if (state_vector.region(pid) == state_vector.face(pid))
        {
            db = check_shell(r, omega,
                             state_vector.face(pid),
                             state_vector.face(pid),
                             state_vector.region(pid) + 1,
                             state_vector.face(pid),
                             state_vector, pid);

            // if we can't hit the shell because of a glancing shot + floating
            // point error, update the region since we won't traverse the
            // shell
            if (db < 0.0)
            {
                ++state_vector.region(pid);
            }
        }

        // now check for hitting shells greater than the current shell as long
        // as we aren't on the last face
        if (state_vector.face(pid) < d_num_shells - 1)
        {
            check_shell(r, omega, state_vector.face(pid) + 1, Geo_State_t::NONE,
                        state_vector.region(pid) + 1,
                        state_vector.face(pid) + 1,
                        state_vector, pid);
        }

        // now check for hitting shells less than the current shell as long as
        // we aren't on the first face
        if (state_vector.face(pid) > 0)
        {
            check_shell(r, omega, state_vector.face(pid) - 1, Geo_State_t::NONE,
                        state_vector.region(pid) - 1,
                        state_vector.face(pid) - 1,
                        state_vector, pid);
        }
    }

    // otherwise, find out where we are and check the bounding shells
    else
    {
        // check for hitting lowest shell
        if (state_vector.region(pid) == 0)
        {
            // we can only hit the lowest shell
            check_shell(r, omega, 0, Geo_State_t::NONE, 1, 0, state_vector, pid);
        }

        // check for hitting highest shell
        else if (state_vector.region(pid) == d_mod_region)
        {
            // we can only hit the outer shell
            check_shell(r, omega,
                        d_num_shells - 1,
                        Geo_State_t::NONE,
                        d_num_shells - 1,
                        d_num_shells - 1,
                        state_vector, pid);
        }

        // otherwise we are between shells
        else
        {
            DEVICE_CHECK(state_vector.region(pid) - 1 >= 0);

            // check hitting lower shell
            check_shell(r, omega,
                        state_vector.region(pid) - 1,
                        Geo_State_t::NONE,
                        state_vector.region(pid) - 1,
                        state_vector.region(pid) - 1,
                        state_vector, pid);

            // check hitting higher shell
            check_shell(r, omega,
                        state_vector.region(pid),
                        Geo_State_t::NONE,
                        state_vector.region(pid) + 1,
                        state_vector.region(pid),
                        state_vector, pid);
        }
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Distance to external radial surfaces.
 */
__device__
void RTK_Cell::dist_to_radial_face(
    int                 axis,
    double              p,
    double              dir,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    double db = profugus::constants::huge;
    int face  = 0;

    // check high/low faces
    if (dir > 0.0)
    {
        db   = (d_extent[axis][HI] - p) / dir;
        face = Geo_State_t::plus_face(axis);
    }
    else if (dir < 0.0)
    {
        db   = (d_extent[axis][LO] - p) / dir;
        face = Geo_State_t::minus_face(axis);
    }
    DEVICE_CHECK(db >= 0.0);

    // updated distance to boundary info
    if (db < state_vector.dist_to_next_region(pid))
    {
        state_vector.dist_to_next_region(pid) = db;
        state_vector.exiting_face(pid)        = face;
        state_vector.next_face(pid)           = Geo_State_t::NONE;
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Distance to external radial surfaces.
 */
__device__
void RTK_Cell::dist_to_axial_face(
    double              p,
    double              dir,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    double db = profugus::constants::huge;
    int face  = 0;

    // check high/low faces
    if (dir > 0.0)
    {
        db   = (d_z - p) / dir;
        face = Geo_State_t::PLUS_Z;
    }
    else if (dir < 0.0)
    {
        db   = -p / dir;
        face = Geo_State_t::MINUS_Z;
    }
    DEVICE_CHECK(db >= 0.0);

    // updated distance to boundary info
    if (db < state_vector.dist_to_next_region(pid))
    {
        state_vector.dist_to_next_region(pid) = db;
        state_vector.exiting_face(pid)        = face;
        state_vector.next_face(pid)           = Geo_State_t::NONE;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Distance to vessel.
 */
__device__
void RTK_Cell::dist_to_vessel(
    const Space_Vector &r,
    const Space_Vector &omega,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    using def::X; using def::Y;

    // return if there is no vessel
    if (!d_vessel)
    {
        return;
    }

    DEVICE_REQUIRE(d_num_shells == 0);

    // local boolean for hitting a vessel wall
    bool hit = false;

    // distance to boundary
    double db = profugus::constants::huge;

    // check distance to first shell
    if (d_inner)
    {
        // only check if we aren't currently on the vessel face
        if (state_vector.face(pid) != Geo_State_t::R0_VESSEL)
        {
            db = dist_to_shell(
                l2g(r[X], X), l2g(r[Y], Y), omega[X], omega[Y], d_R0,
                Geo_State_t::R0_VESSEL);
        }

        // update the distance to boundary
        if (db > 0.0)
        {
            if (db < state_vector.dist_to_next_region(pid))
            {
                state_vector.dist_to_next_region(pid) = db;
                state_vector.next_face(pid)           = Geo_State_t::R0_VESSEL;
                state_vector.exiting_face(pid)        = Geo_State_t::INTERNAL;
                hit                       = true;
            }
        }
    }

    // check distance to second shell
    if (d_outer)
    {
        // only check if we aren't currently on the vessel face
        if (state_vector.face(pid) != Geo_State_t::R1_VESSEL)
        {
            db = dist_to_shell(
                l2g(r[X], X), l2g(r[Y], Y), omega[X], omega[Y], d_R1,
                Geo_State_t::R1_VESSEL);
        }

        // update the distance to boundary
        if (db > 0.0)
        {
            if (db < state_vector.dist_to_next_region(pid))
            {
                state_vector.dist_to_next_region(pid) = db;
                state_vector.next_face(pid)           = Geo_State_t::R1_VESSEL;
                state_vector.exiting_face(pid)        = Geo_State_t::INTERNAL;
                hit                       = true;
            }
        }
    }

    // if we are in the vessel and we hit a surface turn off the vessel;
    // otherwise we enter the vessel
    if (hit)
    {
        if (state_vector.region(pid) == Geo_State_t::VESSEL)
        {
            state_vector.next_region(pid) = Geo_State_t::MODERATOR;
        }
        else
        {
            state_vector.next_region(pid) = Geo_State_t::VESSEL;
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Distance to a shell.
 */
__device__
double RTK_Cell::dist_to_shell(
    double x,
    double y,
    double omega_x,
    double omega_y,
    double r,
    int    face) const
{
    // initialize distance to boundary
    double db = -1.0;

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
            db = d2;
        else if (d2 < 0.0)
            db = d1;
        else if (face < d_num_shells)
            db = max(d1, d2);
        else
            db = min(d1, d2);
    }

    return db;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update state if it hits a shell.
 */
__device__
double RTK_Cell::check_shell(
    const Space_Vector &r,
    const Space_Vector &omega,
    int                 shell,
    int                 face,
    int                 next_region,
    int                 next_face,
    Geo_State_Vector_t &state_vector,
    int                 pid) const
{
    using def::X; using def::Y; using def::Z;

    DEVICE_REQUIRE(shell >= 0 && shell < d_num_shells);

    // calculate the distance to the requested shell
    double db = dist_to_shell(r[X], r[Y], omega[X], omega[Y], d_r[shell], face);

    // check the distance to boundary
    //    a) if it intersects the shell, and
    //    b) if it is the smallest distance
    if (db > 0.0)
    {
        if (db < state_vector.dist_to_next_region(pid))
        {
            state_vector.dist_to_next_region(pid) = db;
            state_vector.next_region(pid)         = next_region;
            state_vector.next_face(pid)           = next_face;
            state_vector.exiting_face(pid)        = Geo_State_t::INTERNAL;
        }
    }

    return db;
}

} // end namespace cuda_profugus

#endif // MC_cuda_rtk_RTK_Cell_i_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Cell.i.cuh
//---------------------------------------------------------------------------//
