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

namespace gpu
{

//===========================================================================//
/*!
 * \class Geometry
 * \brief ACC Geometry wrapper.
 */
//===========================================================================//

class Geometry
{
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
#pragma acc routine seq
    int num_cells() const { return d_N[0]*d_N[1]*d_N[2]; }
};

} // end namespace gpu

#endif // acc_Geometry_hh

//---------------------------------------------------------------------------//
//                 end of Geometry.hh
//---------------------------------------------------------------------------//
