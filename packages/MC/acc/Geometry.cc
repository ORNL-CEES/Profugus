//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   acc/Geometry.cc
 * \author Thomas M. Evans
 * \date   Tue Oct 28 16:37:36 2014
 * \brief  Geometry member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "Geometry.hh"

namespace gpu
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Geometry::Geometry(int    N,
                   double d)
    : d_edges(3, std::vector<double>(N+1, 0.0))
{
    // build the grid edges
    for (int n = 1; n < N+1; ++n)
    {
        d_edges[0][n] = d_edges[0][n-1] + d;
        d_edges[1][n] = d_edges[1][n-1] + d;
        d_edges[2][n] = d_edges[2][n-1] + d;
    }

    d_x = &d_edges[0][0];
    d_y = &d_edges[1][0];
    d_z = &d_edges[2][0];

    std::fill(std::begin(d_N), std::end(d_N), N);

#pragma acc enter data copyin(this)
#pragma acc enter data copyin(d_x[0:N+1], d_y[0:N+1], d_z[0:N+1])
#pragma acc enter data copyin(d_N[0:3])
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Geometry::~Geometry()
{
#pragma acc exit data delete(this)
#pragma acc exit data delete(d_x, d_y, d_z, d_N)
}

} // end namespace gpu

//---------------------------------------------------------------------------//
//                 end of Geometry.cc
//---------------------------------------------------------------------------//
