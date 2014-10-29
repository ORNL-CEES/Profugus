//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   acc/Geometry.hh
 * \author Thomas M. Evans
 * \date   Tue Oct 28 16:37:36 2014
 * \brief  Geometry class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef acc_Geometry_hh
#define acc_Geometry_hh

#include <vector>
#include "core/geometry/Definitions.hh"
#include "core/mc/Definitions.hh"

namespace acc
{

//===========================================================================//
/*!
 * \struct Mesh_State
 * \brief Track particle location on a cartesian mesh
 */
//===========================================================================//

struct Geometry_State
{
  public:
    // >>> PUBLIC DATA

    //! Indices along the mesh grid if inside (invalid if not)
    int ijk[3];

    //! Position
    double pos[3];

    //! Direction
    double dir[3];

    //! Coordinates of next cell (not pickled)
    int next_ijk[3];

    //! Distance to next cell
    double next_dist;
};

//===========================================================================//
/*!
 * \class Geometry
 * \brief ACC Geometry wrapper.
 */
//===========================================================================//

class Geometry
{
  public:
    typedef profugus::events::Event Event_Type;

  private:
    // >>> CPU DATA

    // Edges.
    std::vector<std::vector<double>> d_edges;

    // Matids.
    std::vector<int> d_matids;

    // Boundaries.
    int d_bnds[6];

    // >>> DEVICE DATA

    // Edges.
    double *d_x, *d_y, *d_z;
    int     d_N[3];

    // Matids.
    int *d_m;

    // Boundaries.
    int d_b[6];

    // Work arrays on the GPU.
    double d_work[3];

  public:
    // Constructor.
    Geometry(int N, double d, const std::vector<int> &matids,
             const int *bnds);

    // Destructor.
    ~Geometry();

    //! Get number of cells.
    int num_cells() const { return d_N[0]*d_N[1]*d_N[2]; }

    // Initialize the geometry state.
#pragma acc routine seq
    void initialize(const double *r, const double *direction,
                    Geometry_State &state);

    // Get distance to boundary.
#pragma acc routine seq
    double distance_to_boundary(Geometry_State &state) const;

    // Change the direction through an angle.
#pragma acc routine seq
    void change_direction(double costheta, double phi,
                          Geometry_State& state) const;

    // Reflect the direction at a reflecting surface.
#pragma acc routine seq
    bool reflect(Geometry_State& state) const;

    // Return the boundary state.
#pragma acc routine seq
    profugus::geometry::Boundary_State boundary_state(
        const Geometry_State &state) const;

    //! Move to and cross a surface in the current direction.
    void move_to_surface(Geometry_State& state) const
    {
        state.pos[0] += state.dir[0] * state.next_dist;
        state.pos[1] += state.dir[1] * state.next_dist;
        state.pos[2] += state.dir[2] * state.next_dist;

        state.ijk[0] = state.next_ijk[0];
        state.ijk[1] = state.next_ijk[1];
        state.ijk[2] = state.next_ijk[2];
    }

    //! Move a distance \e d to a point in the current direction.
    void move_to_point(double d, Geometry_State& state) const
    {
        state.pos[0] += state.dir[0] * state.next_dist;
        state.pos[1] += state.dir[1] * state.next_dist;
        state.pos[2] += state.dir[2] * state.next_dist;
    }

    //! Return the current position.
    const double* position(const Geometry_State& state) const
    {
        return state.pos;
    }

    //! Return the current direction.
    const double* direction(const Geometry_State& state) const
    {
        return state.dir;
    }

    //! Return the current material ID
    profugus::geometry::matid_type matid(const Geometry_State& state) const
    {
        return 0;
    }
};

} // end namespace acc

#endif // acc_Geometry_hh

//---------------------------------------------------------------------------//
//                 end of Geometry.hh
//---------------------------------------------------------------------------//
