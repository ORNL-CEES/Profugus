//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_Cell.cuh
 * \author Thomas Evans
 * \date   Mon Nov 21 14:20:32 2016
 * \brief  RTK_Cell kernel declarations.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_RTK_Cell_cuh
#define MC_cuda_rtk_RTK_Cell_cuh

#include <thrust/device_vector.h>

#include "Utils/utils/Constants.hh"
#include "CudaUtils/cuda_utils/CudaDBC.hh"
#include "CudaUtils/cuda_utils/Device_View_Field.hh"
#include "CudaUtils/cuda_utils/Device_Vector_Lite.hh"
#include "CudaUtils/cuda_utils/Device_Memory_Manager.hh"
#include "CudaUtils/cuda_utils/Utility_Functions.hh"
#include "MC/geometry/RTK_Cell.hh"
#include "RTK_State.cuh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class RTK_Cell
 * \brief CUDA implementation of profugus::RTK_Cell.
 */
//===========================================================================//

class RTK_Cell
{
  public:
    //@{
    //! Types.
    using Geo_State_t  = RTK_State;
    using Space_Vector = Geo_State_t::Space_Vector;
    //@}

  private:
    // Types
    using View_Int = cuda::const_Device_View_Field<int>;
    using View_Dbl = cuda::const_Device_View_Field<double>;

  private:
    // >>> DATA

    // Moderator id.
    int d_mod_id;

    // Shell radii and material ids.
    View_Dbl d_r;
    View_Int d_ids;

    // Radial dimensions (pitch).
    cuda_utils::Device_Vector_Lite<double, 2> d_xy;

    // Z-extent.
    double d_z;

  public:
    // Constructor.
    RTK_Cell(int mod_id, View_Dbl r, View_Int ids, double (&extents)[2][2],
             double height, int num_segments);

    // Add a vessel.
    void add_vessel(int vessel_id, double R0, double R1, double (&offsets)[2]);

    // Initialize a state.
    __device__
    void initialize(const Space_Vector &r, Geo_State_t &state) const;

    // Track to next boundary.
    __device__
    void distance_to_boundary(const Space_Vector &r,
                              const Space_Vector &omega,
                              Geo_State_t &state);

    // Update a state at collision sites.
    __device__
    void update_state(Geo_State_t &state) const;

    // Cross a surface.
    __device__
    void cross_surface(Geo_State_t &state) const;

    // Query to find region.
    __device__
    int region(double x, double y) const;

    // Query to find segment.
    __device__
    int segment(double x, double y) const;

    // Query to find cell.
    __device__
    inline int cell(int region, int segment) const;

    // Get extents of this geometry element in the parent reference frame
    __device__
    inline void get_extents(Space_Vector &lower, Space_Vector &upper) const;

    // Return the material id for a region in the pin-cell.
    __device__
    inline int matid(int region) const;

    //! Return the number of regions.
    __device__
    int num_regions() const { return d_num_regions; }

    //! Return the number of shells.
    __device__
    int num_shells() const { return d_num_shells; }

    //! Return the number of segments.
    __device__
    int num_segments() const { return d_segments; }

    //! Number of cells.
    __device__
    int num_cells() const { return d_num_cells; }

    //! Return pin-pitch in X or Y dimension.
    __device__
    double pitch(int d) const { DEVICE_REQUIRE(d < 2); return d_xy[d]; }

    //! Return height.
    __device__
    double height() const { return d_z; }

    //! Return shell radii.
    __device__
    View_Dbl radii() const { return d_r; }

    //! Query for vessel.
    __device__
    bool has_vessel() const { return d_vessel; }

  private:
    // >>> IMPLEMENTATION

    // Intersections with shells.
    __device__
    void calc_shell_db(const Space_Vector &r, const Space_Vector &omega,
                       Geo_State_t &state);

    // Distance to external surface.
    __device__
    void dist_to_radial_face(int axis, double p, double dir,
                             Geo_State_t &state);
    __device__
    void dist_to_axial_face(double p, double dir, Geo_State_t &state);

    // Distance to vessel.
    __device__
    void dist_to_vessel(const Space_Vector &r, const Space_Vector &omega,
                        Geo_State_t &state);

    // Distance to a shell.
    __device__
    double dist_to_shell(double x, double y, double omega_x, double omega_y,
                         double r, int face);

    // Update state if it hits a shell.
    __device__
    double check_shell(const Space_Vector &r, const Space_Vector &omega,
                       int shell, int face, int next_region, int next_face,
                       Geo_State_t &state);

    // Transform to vessel coordinates.
    __host__ __device__
    double l2g(double local, int dir) const
    {
        return local + d_offsets[dir];
    }

    // Low/High enum.
    enum LOHI { LO = 0, HI = 1 };

    // Half-pitch.
    double d_extent[2][2];

    // Number of regions and shells.
    int d_num_shells;
    int d_num_regions;

    // Number of segments.
    int d_segments;
    int d_seg_faces;

    // Number of internal faces (segment faces + shells).
    int d_num_int_faces;

    // Moderator region id.
    int d_mod_region;

    // Number of cells.
    int d_num_cells;

    // Vessel parameters.
    bool   d_vessel;         // indicates this cell has a vessel
    double d_offsets[2];     // radial offsets from origin of outer rtk-array
                             // to origin of the pincell
    double d_R0, d_R1;       // inner and outer vessel radii
    bool   d_inner, d_outer; // booleans for inner/outer vessel bisections
    int    d_vessel_id;      // vessel mat id
};

//===========================================================================//
/*!
 * \class RTK_Cell_DMM
 * \brief Device_Memory_Manager for the device-based RTK_Cell class.
 */
//===========================================================================//

class RTK_Cell_DMM
    : public cuda::Device_Memory_Manager<RTK_Cell>
{
    using Base = cuda::Device_Memory_Manager<RTK_Cell>;

  public:
    //! Host type.
    using Host_RTK_Cell = profugus::RTK_Cell;

  private:
    // >>> DEVICE MEMORY MANAGEMENT DATA

    // Shell radii.
    thrust::device_vector<double> d_r;

    // Shell ids.
    thrust::device_vector<int> d_ids;

  public:
    // Constructor.
    RTK_Cell_DMM(const Host_RTK_Cell &host_cell);

    // >>> DERIVED INTERFACE

    // Create a device instance.
    RTK_Cell device_instance();

  private:
    // >>> IMPLEMENTATION

    // POD data need to construct device object.
    int    d_mod_id;
    int    d_num_segments;
    double d_z;
    double d_extent[2][2];

    // Vessel data
    bool   d_vessel;
    int    d_vessel_id;
    double d_offsets[2];
    double d_R0, d_R1;
};

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// CUDA FUNCTIONS
//---------------------------------------------------------------------------//

#include "RTK_Cell.i.cuh"

#endif // MC_cuda_rtk_RTK_Cell_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Cell.cuh
//---------------------------------------------------------------------------//
