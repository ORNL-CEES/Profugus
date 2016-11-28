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
// PUBLIC FUNCTIONS
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
/*
 * \brief Distance to external radial surfaces.
 */
__device__
void RTK_Cell::dist_to_radial_face(
    int          axis,
    double       p,
    double       dir,
    Geo_State_t &state)
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
    if (db < state.dist_to_next_region)
    {
        state.dist_to_next_region = db;
        state.exiting_face        = face;
        state.next_face           = Geo_State_t::NONE;
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Distance to external radial surfaces.
 */
__device__
void RTK_Cell::dist_to_axial_face(
    double       p,
    double       dir,
    Geo_State_t &state)
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
    if (db < state.dist_to_next_region)
    {
        state.dist_to_next_region = db;
        state.exiting_face        = face;
        state.next_face           = Geo_State_t::NONE;
    }
}

} // end namespace cuda_profugus

#endif // MC_cuda_rtk_RTK_Cell_i_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Cell.i.cuh
//---------------------------------------------------------------------------//
