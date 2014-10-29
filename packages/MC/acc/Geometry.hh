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

    //! Faces in pin-cell.
    enum Face {
        MINUS_X = 0,
        PLUS_X ,
        MINUS_Y,
        PLUS_Y ,
        MINUS_Z,
        PLUS_Z,
        END_FACES
    };

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

    // >>> DEVICE DATA

    // Edges.
    double *d_x, *d_y, *d_z;
    int     d_N[3];

  public:
    // Constructor.
    Geometry(int N, double d);

    // Destructor.
    ~Geometry();

    //! Get number of cells.
    int num_cells() const { return d_N[0]*d_N[1]*d_N[2]; }

    // Get distance to boundary.
#pragma acc routine seq
    double distance_to_boundary(Geometry_State &state) const;
};

} // end namespace acc

#endif // acc_Geometry_hh

//---------------------------------------------------------------------------//
//                 end of Geometry.hh
//---------------------------------------------------------------------------//
