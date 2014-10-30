//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   acc/Geometry_State.hh
 * \author Seth R Johnson
 * \date   Thu Oct 30 13:27:28 2014
 * \brief  Geometry_State class declaration.
 * \note   Copyright (c) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef acc_Geometry_State_hh
#define acc_Geometry_State_hh

namespace acc
{

//===========================================================================//
/*!
 * \struct Geometry_State
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

//---------------------------------------------------------------------------//
} // end namespace acc

#endif // acc_Geometry_State_hh

//---------------------------------------------------------------------------//
// end of acc/Geometry_State.hh
//---------------------------------------------------------------------------//
