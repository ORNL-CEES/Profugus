//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/RTK_Cell.hh
 * \author Thomas M. Evans
 * \date   Tuesday April 29 15:40:59 2014
 * \brief  RTK_Cell class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_RTK_Cell_hh
#define geometry_RTK_Cell_hh

#include <ostream>
#include <vector>

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "utils/Vector_Lite.hh"
#include "RTK_State.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class RTK_Cell
 * \brief Defines a single pin-cell geometry for MC transport.
 *
 * The following figure represents a pin-cell:
 * \image html pin_cell.png "Pin cell"
 * This pin cells has: \e 4-shells; \e 5-regions; \e 4-segments; and \e
 * 20-cells. The segments are indicated by the red-dashed lines.  The cell
 * indexing is:
 * \f[
   \mbox{cell\_index} = \mbox{region} + \mbox{segment}\times N_{r}
 * \f]
 * where \f$N_r\f$ is the number of regions.  The pin-cell construction is
 * constrained such that the number of regions per segment is always a
 * constant.  With these definitions, the number of cells in each pin-cell is:
 * \f[
   N_{c} = N_{r}\times N_{s}
 * \f]
 *
 * You can also specify gap on up to 2 sides of the pin, thus representing
 * assembly gap.  There are no explicit "gap" cells, the gap is just included
 * as extended moderator region on the appropriate side.  For example,
 * specifying a pin with pitch 1.26 with a gap on the right side of 0.1 would
 * return an "effective" pitch of 1.36 in the x-direction.  The pins are still
 * centered about x = 0.0:
 * \verbatim
    -0.63       0.0         0.73
      |          .           |
  left-edge    origin   right-edge
   \endverbatim
 *
 * The relationship between internal faces, regions, and cells is illustrated
 * in the following figure:
 * \image html pin_cell_test.png "Pin-cell faces, regions, and * cells."
 * The regions are indicated with blue ids, the faces are indicated with
 * red ids, and the cells each are given a separate color.  In this example
 * there are
 * - 4 regions
 * - 3 shells (with 3 faces)
 * - 4 segments (requiring 2 segment faces)
 * - 5 total internal faces (3 shell faces + 2 segment faces)
 * - 16 cells (4 segments times 4 regions)
 * .
 * See geometry/rtk/test/tstRTK_Cell.cc for a tracking example through this pin
 * cell configuration.
 *
 * There are 2 types of "special" pins: homogeneous and vessel.  Homogeneous
 * pins are simple boxes.  Vessel pins are homogeneous cells with a core
 * vessel bisecting it as shown below.
 * \image html core_vessel.png "Homogeneous and vessel pins."
 */
/*!
 * \example geometry/test/tstRTK_Cell.cc
 *
 * Test of RTK_Cell.
 */
//===========================================================================//

class RTK_Cell
{
  public:
    //@{
    //! Useful typedefs.
    typedef def::Vec_Int           Vec_Int;
    typedef def::Vec_Dbl           Vec_Dbl;
    typedef def::Space_Vector      Space_Vector;
    typedef Vector_Lite<double, 4> Gap_Vector;
    typedef RTK_State              Geo_State_t;
    //@}

  private:
    // >>> DATA

    // Moderator id.
    int d_mod_id;

    // Shell radii and material ids.
    Vec_Dbl d_r;
    Vec_Int d_ids;

    // Radial dimensions (pitch).
    Vector_Lite<double, 2> d_xy;

    // Z-extent.
    double d_z;

  public:
    // Constructors.
    RTK_Cell(int mod_id, double pitch, double height, int segments = 1);
    RTK_Cell(int mod_id, double dx, double dy, double height, int segments = 1);
    RTK_Cell(int fuel_id, double r, int mod_id, double pitch, double height,
             int segments = 1);
    RTK_Cell(const Vec_Int &ids, const Vec_Dbl &r, int mod_id, double pitch,
             double height, int segments = 1);
    RTK_Cell(const Vec_Int &ids, const Vec_Dbl &r, int mod_id, double pitch,
             double height, const Gap_Vector &gap, int segments = 1);
    RTK_Cell(int mod_id, double pitch, double height, const Gap_Vector &gap,
             int segments = 1);
    RTK_Cell(int mod_id, double dx, double dy, double height, double R0,
             double R1, double x_offset, double y_offset, int vessel_id);

    // Pin-cells are completed on construction.
    bool completed() const { return true; }

    // Add vessel to pin-cell.
    bool add_vessel(double R0, double R1, double x_offset, double y_offset,
                    int vessel_id);

    // Initialize a state.
    void initialize(const Space_Vector &r, Geo_State_t &state) const;

    // Track to next boundary.
    void distance_to_boundary(const Space_Vector &r,
                              const Space_Vector &omega,
                              Geo_State_t &state);

    // Update a state at collision sites.
    void update_state(Geo_State_t &state) const;

    // Cross a surface.
    void cross_surface(Geo_State_t &state) const;

    // Query to find region.
    int region(double x, double y) const;

    // Query to find segment.
    inline int segment(double x, double y) const;

    // Query to find cell.
    inline int cell(int region, int segment) const;

    // Get extents of this geometry element in the parent reference frame
    inline void get_extents(Space_Vector &lower, Space_Vector &upper) const;

    // Return the material id for a region in the pin-cell.
    inline int matid(int region) const;

    //! Return the number of regions.
    int num_regions() const { return d_num_regions; }

    //! Return the number of shells.
    int num_shells() const { return d_num_shells; }

    //! Return the number of segments.
    int num_segments() const { return d_segments; }

    //! Number of cells.
    int num_cells() const { return d_num_cells; }

    //! Return pin-pitch in X or Y dimension.
    double pitch(int d) const { REQUIRE(d < 2); return d_xy[d]; }

    //! Return height.
    double height() const { return d_z; }

    //! Return shell radii.
    const Vec_Dbl& radii() const { return d_r; }

    // Output pin for diagnostics.
    void output(std::ostream &ostream, int array_id, int pin_id) const;

    // Query for vessel.
    bool has_vessel() const { return d_vessel; }

    // Vessel radius and offsets.
    bool vessel_data(double &R0, double &R1, double &xc, double &yc) const;

  private:
    // >>> IMPLEMENTATION

    // Intersections with shells.
    void calc_shell_db(const Space_Vector &r, const Space_Vector &omega,
                       Geo_State_t &state);

    // Distance to external surface.
    inline void dist_to_radial_face(int axis, double p, double dir,
                                    Geo_State_t &state);
    inline void dist_to_axial_face(double p, double dir, Geo_State_t &state);

    // Distance to vessel.
    void dist_to_vessel(const Space_Vector &r, const Space_Vector &omega,
                        Geo_State_t &state);

    // Distance to a shell.
    void dist_to_shell(double x, double y, double omega_x, double omega_y,
                       double r, int face);

    // Update state if it hits a shell.
    void check_shell(const Space_Vector &r, const Space_Vector &omega,
                     int shell, int face, int next_region, int next_face,
                     Geo_State_t &state);

    // Transform to vessel coordinates.
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

    // Work variables.
    double d_db;
    int    d_face;
    int    d_segment;

    // Vessel parameters.
    bool d_vessel;         // indicates this cell has a vessel
    double d_offsets[2];   // radial offsets from origin of outer rtk-array to
                           // origin of the pincell
    double d_R0, d_R1;     // inner and outer vessel radii
    bool d_inner, d_outer; // booleans for inner/outer vessel bisections
    int d_vessel_id;       // vessel mat id
};

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "RTK_Cell.i.hh"

#endif // geometry_RTK_Cell_hh

//---------------------------------------------------------------------------//
//              end of geometry/RTK_Cell.hh
//---------------------------------------------------------------------------//
