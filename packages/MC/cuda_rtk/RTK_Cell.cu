//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_Cell.cu
 * \author Tom Evans
 * \date   Mon Nov 28 12:33:05 2016
 * \brief  RTK_Cell kernel definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <cmath>

#include "RTK_Cell.cuh"

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
// RTK_CELL MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
RTK_Cell::RTK_Cell(int      mod_id,
                   View_Dbl r,
                   View_Int ids,
                   double   (&extents)[2][2],
                   double   height,
                   int      num_segments)
    : d_mod_id(mod_id)
    , d_r(r)
    , d_ids(ids)
    , d_z(height)
    , d_num_shells(r.size())
    , d_num_regions(d_num_shells + 1)
    , d_segments(num_segments)
    , d_seg_faces(d_segments / 2)
    , d_num_int_faces(d_seg_faces + d_num_shells)
    , d_mod_region(d_num_regions - 1)
    , d_num_cells(d_num_regions * d_segments)
    , d_vessel(false)
{
    using def::X; using def::Y;

    // Set the extents
    d_extent[X][LO] = extents[X][LO];
    d_extent[X][HI] = extents[X][HI];
    d_extent[Y][LO] = extents[Y][LO];
    d_extent[Y][HI] = extents[Y][HI];

    // Set the pitch
    d_xy[X] = d_extent[X][HI] - d_extent[X][LO];
    d_xy[Y] = d_extent[Y][HI] - d_extent[Y][LO];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add a vessel.
 */
void RTK_Cell::add_vessel(
    int    vessel_id,
    double R0,
    double R1,
    double (&offsets)[2])
{
    using def::X; using def::Y;

    // Set internal data
    d_vessel     = true;
    d_vessel_id  = vessel_id;
    d_offsets[0] = offsets[0];
    d_offsets[1] = offsets[1];
    d_R0         = -1.0;
    d_R1         = -1.0;
    d_inner      = false;
    d_outer      = false;

    // Near and far corners relative to the origin of vessel cylinder
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
// RTK_CELL_DMM MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
RTK_Cell_DMM::RTK_Cell_DMM(
    const Host_RTK_Cell &host_cell)
    : d_num_segments(host_cell.num_segments())
    , d_z(host_cell.height())
{
    using def::X; using def::Y;

    // Number of shells in this cell
    auto num_shells = host_cell.num_shells();

    // Calculate radii and ids
    if (num_shells > 0)
    {
        // Radii
        const auto &radii = host_cell.radii();
        CHECK(!radii.empty());

        d_r = thrust::device_vector<double>(radii.begin(), radii.end());

        // IDS
        std::vector<int> ids(num_shells);
        for (int n = 0; n < num_shells; ++n)
        {
            ids[n] = host_cell.matid(n);
        }
        d_ids = thrust::device_vector<int>(ids.begin(), ids.end());
    }

    // Number of regions in this cell
    auto num_regions = host_cell.num_regions();
    CHECK(num_regions == num_shells + 1);

    // Store the moderator id
    d_mod_id = host_cell.matid(num_regions - 1);

    // Get the extents
    Host_RTK_Cell::Space_Vector low, high;
    host_cell.get_extents(low, high);
    d_extent[X][0] = low[X];
    d_extent[X][1] = high[X];
    d_extent[Y][0] = low[Y];
    d_extent[Y][1] = high[Y];
    CHECK(high[2] == d_z);

    // Store vessel data
    d_vessel = host_cell.vessel_data(
        d_vessel_id, d_R0, d_R1, d_offsets[0], d_offsets[1]);

    ENSURE(host_cell.num_cells() == num_regions * d_num_segments);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct a device RTK_Cell object.
 */
RTK_Cell RTK_Cell_DMM::device_instance()
{
    RTK_Cell cell(d_mod_id,
                  cuda::make_view(d_r),
                  cuda::make_view(d_ids),
                  d_extent,
                  d_z,
                  d_num_segments);

    // Add the vessel if it exists
    if (d_vessel) cell.add_vessel(d_vessel_id, d_R0, d_R1, d_offsets);

    return cell;
}

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Cell.cu
//---------------------------------------------------------------------------//
