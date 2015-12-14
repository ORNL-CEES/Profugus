//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_geometry/Cartesian_Mesh.hh
 * \author Thomas M. Evans
 * \date   Mon Jul 21 16:55:00 2014
 * \brief  Cartesian_Mesh class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_geometry_Cartesian_Mesh_hh
#define cuda_geometry_Cartesian_Mesh_hh

#include <math_constants.h>

#include <vector>
#include <memory>

#include "harness/DBC.hh"
#include "cuda_utils/Definitions.hh"
#include "cuda_utils/Device_Vector.hh"
#include "utils/Definitions.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Cartesian_Mesh
 * \brief Common functionality for Cartesian meshes.
 *
 * Provides mesh-based functionality for tallies.
 */
/*!
 * \example geometry/test/tstCartesian_Mesh.cc
 *
 * Test of Cartesian_Mesh.
 */
//===========================================================================//

class Cartesian_Mesh
{
  public:
    //@{
    //! Mesh typedefs
    typedef int    dim_type;    //!< indexing along a single dimension
    typedef size_t size_type;   //!< indexing for cells

    typedef std::vector<double>                Vec_Dbl;
    typedef cuda::arch::Device                 Arch;
    typedef cuda::Device_Vector<Arch,double>   Dbl_Vec;
    typedef cuda::Device_Vector<Arch,int>      Int_Vec;
    typedef std::shared_ptr<Dbl_Vec>           SP_Dbl_Vec;
    typedef int                                Dim_Vec[3];
    //@}

  private:
    // >>> DATA

    // Cell edges.
    Dbl_Vec d_x_edges_vec;
    Dbl_Vec d_y_edges_vec;
    Dbl_Vec d_z_edges_vec;
    double *d_x_edges;
    double *d_y_edges;
    double *d_z_edges;
    int d_cells_x;
    int d_cells_y;
    int d_cells_z;

    // Number of global cells
    size_type d_num_cells;

    // Dimensionality
    dim_type d_dimension;

    SP_Dbl_Vec d_volumes_vec;
    double *d_volumes;

  public:
    // Construct from xyz edges.
    Cartesian_Mesh(const Vec_Dbl& x_edges, const Vec_Dbl& y_edges,
                   const Vec_Dbl& z_edges);

    // >>> ACCESSORS

    //! Get number of cells.
    __host__ __device__ size_type num_cells() const { return d_num_cells; }

    //! Number of cells along an axis
    __device__ dim_type num_cells_along(dim_type d) const
    {
        if( d == def::I )
            return d_cells_x;
        else if( d == def::J )
            return d_cells_y;
        else if( d == def::K )
            return d_cells_z;
        return -1;
    }

    //! Dimension of mesh.
    __host__ __device__ dim_type dimension() const { return d_dimension; }

    //! Return cell edges along a given direction.
    __device__ double * edges(dim_type d) const
    {
        if( d == def::I )
            return d_x_edges;
        else if( d == def::J )
            return d_y_edges;
        else if( d == def::K )
            return d_z_edges;
        return 0;
    }

    // >>> INDEX CONVERSION

    // Convert cardinal index to (i,j) or (i,j,k).
    __device__ void cardinal(
        size_type cell, dim_type& i, dim_type& j, dim_type& k) const
    {
        k = cell / (d_cells_x * d_cells_y);
        cell = cell % (d_cells_x * d_cells_y);
        j = cell / d_cells_x;
        i = cell % d_cells_x;
    }

    // Convert (i,j,k) to cell index
    __host__ __device__ inline bool index(dim_type   i,
                                 dim_type   j,
                                 dim_type   k,
                                 size_type &cell) const
    {
        if( i < d_cells_x && j < d_cells_y && k < d_cells_z )
        {
            cell = i + d_cells_x * (j + d_cells_y * k);
            return true;
        }
        cell = static_cast<size_type>(-1);
        return false;
    }

    // >>> VOLUME CALCULATION

    //! Calculate volume from the global cell id
    __device__ inline double volume(size_type global_cell_index) const
    {
        return d_volumes[global_cell_index];
    }

    // >>> SPATIAL LOCATION

    // Locate a coordinate along a single axis
    __device__ inline dim_type find_upper(
        double r, dim_type axis) const;

    //! Low corner of mesh in \e (i,j,k) direction.
    __device__ double low_corner(dim_type d) const
    {
        if( d == def::I )
            return d_x_edges[0];
        else if( d == def::J )
            return d_y_edges[0];
        else if( d == def::K )
            return d_z_edges[0];
        return CUDART_NAN;
    }

    //! High corner of mesh in \e (i,j,k) direction.
    __device__ double high_corner(dim_type d) const
    {
        if( d == def::I )
            return d_x_edges[d_cells_x];
        else if( d == def::J )
            return d_y_edges[d_cells_y];
        else if( d == def::K )
            return d_z_edges[d_cells_z];
        return CUDART_NAN;
    }
};

} // end namespace cuda_profugus


#endif // cuda_geometry_Cartesian_Mesh_hh

//---------------------------------------------------------------------------//
//                 end of Cartesian_Mesh.hh
//---------------------------------------------------------------------------//
