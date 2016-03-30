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

#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/CudaMacros.hh"
#include "cuda_utils/Definitions.hh"
#include "cuda_utils/Device_Vector.hh"
#include "cuda_utils/Utility_Functions.hh"
#include "cuda_utils/Host_Vector.hh"
#include "geometry/Definitions.hh"
#include "utils/Definitions.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Cartesian_Mesh
 * \brief Common functionality for Cartesian meshes.
 *
 * Provides mesh-based functionality.
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
    typedef int                           dim_type;
    typedef size_t                        size_type;
    typedef profugus::geometry::cell_type cell_type;
    typedef cuda::Space_Vector            Space_Vector;
    typedef cuda::Coordinates             Coordinates;
    typedef std::vector<double>           Vec_Dbl;

    template <class T>
    using Device_Vector = cuda::Device_Vector<cuda::arch::Device,T>;

    template <class T>
    using Host_Vector = cuda::Host_Vector<T>;

    template <class T>
    using SP = std::shared_ptr<T>;
    //@}

  private:
    // >>> DATA

    // Cell edges.
    Device_Vector<double> d_x_edges_vec;
    Device_Vector<double> d_y_edges_vec;
    Device_Vector<double> d_z_edges_vec;

    // On-device pointers
    double *dd_x_edges;
    double *dd_y_edges;
    double *dd_z_edges;

    // Number of cells along each dimension
    int d_cells_x;
    int d_cells_y;
    int d_cells_z;

    // Total number of cells
    cell_type d_num_cells;

    // Dimensionality (always 3 for now)
    dim_type d_dimension;

    // Cell volumes
    SP<Device_Vector<double> > d_volumes_vec;
    SP<Host_Vector<double> > d_volumes_vec_host;

    // On-device pointer
    double *dd_volumes;

  public:

    // Construct from xyz edges.
    Cartesian_Mesh(const Vec_Dbl& x_edges, const Vec_Dbl& y_edges,
                   const Vec_Dbl& z_edges);

    // >>> ACCESSORS

    //! Get number of cells.
    PROFUGUS_HOST_DEVICE_FUNCTION 
    cell_type num_cells() const { return d_num_cells; }

    //! Number of cells along an axis
    PROFUGUS_HOST_DEVICE_FUNCTION 
    dim_type num_cells_along(dim_type d) const
    {
        REQUIRE( d>=0 && d<=3 );
        if( d == def::I )
            return d_cells_x;
        else if( d == def::J )
            return d_cells_y;
        else if( d == def::K )
            return d_cells_z;
        return -1;
    }

    //! Dimension of mesh.
    PROFUGUS_HOST_DEVICE_FUNCTION 
    dim_type dimension() const { return d_dimension; }

    //! Return cell edges along a given direction.
    PROFUGUS_DEVICE_FUNCTION 
    double * edges(dim_type d) const
    {
        REQUIRE( d>=0 && d<=3 );
        if( d == def::I )
            return dd_x_edges;
        else if( d == def::J )
            return dd_y_edges;
        else if( d == def::K )
            return dd_z_edges;
        return 0;
    }

    // >>> INDEX CONVERSION

    // Convert cardinal index to (i,j) or (i,j,k).
    PROFUGUS_HOST_DEVICE_FUNCTION 
    void cardinal(
        cell_type cell, dim_type& i, dim_type& j, dim_type& k) const
    {
        REQUIRE( cell < d_num_cells );
        k = cell / (d_cells_x * d_cells_y);
        ENSURE( k < d_cells_z );
        cell = cell % (d_cells_x * d_cells_y);
        j = cell / d_cells_x;
        ENSURE( j < d_cells_y );
        i = cell % d_cells_x;
        ENSURE( i < d_cells_x );
    }

    // Convert (i,j,k) to cell index
    PROFUGUS_HOST_DEVICE_FUNCTION 
    bool index(dim_type   i,
	       dim_type   j,
	       dim_type   k,
	       cell_type &cell) const
    {
        REQUIRE( i < d_cells_x );
        REQUIRE( j < d_cells_y );
        REQUIRE( k < d_cells_z );
        if( i < d_cells_x && j < d_cells_y && k < d_cells_z )
        {
            cell = i + d_cells_x * (j + d_cells_y * k);
            ENSURE( cell < d_num_cells );
            return true;
        }
        cell = static_cast<cell_type>(-1);
        return false;
    }

    // >>> VOLUME CALCULATION

    //! Calculate volume from the global cell id on the device.
    PROFUGUS_DEVICE_FUNCTION 
    double volume(cell_type global_cell ) const
    {
        REQUIRE( global_cell < d_num_cells );
        return dd_volumes[global_cell];
    }

    //! Calculate volume from the global cell id on the host.
    double volume_host(cell_type global_cell ) const
    {
        REQUIRE( global_cell < d_num_cells );
        return (*d_volumes_vec_host)[global_cell];
    }

    // >>> SPATIAL LOCATION

    // Locate the positon's ijk coordinates with upper edges begin "inside"
    PROFUGUS_DEVICE_FUNCTION 
    void find_upper(const Space_Vector &r, Coordinates &ijk ) const
    {
        ijk.i = find_upper(r.x, def::I);
        ijk.j = find_upper(r.y, def::J);
        ijk.k = find_upper(r.z, def::K);
    }

    // Locate a coordinate along a single axis
    PROFUGUS_DEVICE_FUNCTION 
    dim_type find_upper(
        double r, dim_type axis) const
    {
        REQUIRE( 0 <= axis && axis < d_dimension );
        const double *edges_start;
        const double *edges_end;
        if( axis == def::I )
        {
            edges_start = dd_x_edges;
            edges_end   = dd_x_edges+d_cells_x;
        }
        if( axis == def::J )
        {
            edges_start = dd_y_edges;
            edges_end   = dd_y_edges+d_cells_y;
        }
        if( axis == def::K )
        {
            edges_start = dd_z_edges;
            edges_end   = dd_z_edges+d_cells_z;
        }
        return cuda::utility::lower_bound(edges_start,edges_end,r)
            - edges_start - 1;
    }

    //! Low corner of mesh in \e (i,j,k) direction.
    PROFUGUS_DEVICE_FUNCTION 
    double low_corner(dim_type d) const
    {
        REQUIRE( d>=0 && d<=3 );
	REQUIRE( d == def::I || d == def::J || d == def::K );
        if( d == def::I )
            return dd_x_edges[0];
        else if( d == def::J )
            return dd_y_edges[0];
        else
            return dd_z_edges[0];
    }

    //! High corner of mesh in \e (i,j,k) direction.
    PROFUGUS_DEVICE_FUNCTION 
    double high_corner(dim_type d) const
    {
        REQUIRE( d>=0 && d<=3 );
	REQUIRE( d == def::I || d == def::J || d == def::K );
        if( d == def::I )
            return dd_x_edges[d_cells_x];
        else if( d == def::J )
            return dd_y_edges[d_cells_y];
        else
            return dd_z_edges[d_cells_z];
    }
};

} // end namespace cuda_profugus


#endif // cuda_geometry_Cartesian_Mesh_hh

//---------------------------------------------------------------------------//
//                 end of Cartesian_Mesh.hh
//---------------------------------------------------------------------------//
