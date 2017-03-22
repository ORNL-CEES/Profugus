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
#include <thrust/device_vector.h>

#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/CudaMacros.hh"
#include "cuda_utils/Definitions.hh"
#include "cuda_utils/Utility_Functions.hh"
#include "cuda_utils/Host_Vector.hh"
#include "cuda_utils/Device_Memory_Manager.hh"
#include "cuda_utils/Device_View_Field.hh"
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
    typedef int                                   dim_type;
    typedef size_t                                size_type;
    typedef profugus::geometry::cell_type         cell_type;
    typedef cuda_utils::Space_Vector              Space_Vector;
    typedef cuda_utils::Coordinates               Coordinates;
    typedef cuda::const_Device_View_Field<double> Double_View;
    //@}

  private:
    // >>> DATA

    // Device views
    Double_View d_x_edges;
    Double_View d_y_edges;
    Double_View d_z_edges;

    // Number cells in each dimension
    int d_cells_x;
    int d_cells_y;
    int d_cells_z;

    // Total number of cells
    cell_type d_num_cells;

    // Dimensionality (always 3 for now)
    dim_type d_dimension;

  public:

    // Construct from xyz edges.
    Cartesian_Mesh(Double_View x_edges,
                   Double_View y_edges,
                   Double_View z_edges)
      : d_x_edges(x_edges)
      , d_y_edges(y_edges)
      , d_z_edges(z_edges)
      , d_dimension(3)
    {
        d_cells_x = d_x_edges.size()-1;
        d_cells_y = d_y_edges.size()-1;
        d_cells_z = d_z_edges.size()-1;
        d_num_cells = d_cells_x * d_cells_y * d_cells_z;
    }

    // >>> ACCESSORS

    //! Get number of cells.
    PROFUGUS_HOST_DEVICE_FUNCTION
    cell_type num_cells() const { return d_num_cells; }

    //! Number of cells along an axis
    PROFUGUS_HOST_DEVICE_FUNCTION
    dim_type num_cells_along(dim_type d) const
    {
        DEVICE_REQUIRE( d>=0 && d<=3 );
        if( d == def::I )
            return d_x_edges.size()-1;
        else if( d == def::J )
            return d_y_edges.size()-1;
        else if( d == def::K )
            return d_z_edges.size()-1;
        return -1;
    }

    //! Dimension of mesh.
    PROFUGUS_HOST_DEVICE_FUNCTION
    dim_type dimension() const { return d_dimension; }

    //! Return cell edges along a given direction.
    PROFUGUS_DEVICE_FUNCTION
    Double_View edges(dim_type d) const
    {
        DEVICE_REQUIRE( d>=0 && d<=3 );
        if( d == def::I )
            return d_x_edges;
        else if( d == def::J )
            return d_y_edges;
        else if( d == def::K )
            return d_z_edges;
        return Double_View();
    }

    // >>> INDEX CONVERSION

    // Convert cardinal index to (i,j) or (i,j,k).
    PROFUGUS_HOST_DEVICE_FUNCTION
    void cardinal(
        cell_type cell, dim_type& i, dim_type& j, dim_type& k) const
    {
        DEVICE_REQUIRE( cell < d_num_cells );
        k = cell / (d_cells_x * d_cells_y);
        DEVICE_ENSURE( k < d_cells_z );
        cell = cell % (d_cells_x * d_cells_y);
        j = cell / d_cells_x;
        DEVICE_ENSURE( j < d_cells_y );
        i = cell % d_cells_x;
        DEVICE_ENSURE( i < d_cells_x );
    }

    // Convert (i,j,k) to cell index
    PROFUGUS_HOST_DEVICE_FUNCTION
    bool index(dim_type   i,
	       dim_type   j,
	       dim_type   k,
	       cell_type &cell) const
    {
        DEVICE_REQUIRE( i < d_cells_x );
        DEVICE_REQUIRE( j < d_cells_y );
        DEVICE_REQUIRE( k < d_cells_z );
        if( i < d_cells_x && j < d_cells_y && k < d_cells_z )
        {
            cell = i + d_cells_x * (j + d_cells_y * k);
            DEVICE_ENSURE( cell < d_num_cells );
            return true;
        }
        cell = static_cast<cell_type>(-1);
        return false;
    }

    // >>> SPATIAL LOCATION

    // Locate the positon's ijk coordinates with upper edges begin "inside"
    PROFUGUS_DEVICE_FUNCTION
    void find_upper(const Space_Vector &r, Coordinates &ijk ) const
    {
        using def::I; using def::J; using def::K;
        ijk[I] = find_upper(r[I], I);
        ijk[J] = find_upper(r[J], J);
        ijk[K] = find_upper(r[K], K);
    }

    // Locate a coordinate along a single axis
    PROFUGUS_DEVICE_FUNCTION
    dim_type find_upper(
        double r, dim_type axis) const
    {
        DEVICE_REQUIRE( 0 <= axis && axis < d_dimension );
        const double *edges_start;
        const double *edges_end;
        if( axis == def::I )
        {
            edges_start = d_x_edges.begin();
            edges_end   = d_x_edges.end();
        }
        if( axis == def::J )
        {
            edges_start = d_y_edges.begin();
            edges_end   = d_y_edges.end();
        }
        if( axis == def::K )
        {
            edges_start = d_z_edges.begin();
            edges_end   = d_z_edges.end();
        }
        return cuda_utils::utility::lower_bound(edges_start,edges_end,r)
            - edges_start - 1;
    }
};

//===========================================================================//
/*!
 * \class Cartesian_Mesh_DMM
 * \brief Device memory manager for Cartesian_Mesh
 */
//===========================================================================//

class Cartesian_Mesh_DMM : public cuda::Device_Memory_Manager<Cartesian_Mesh>
{
  public:
    //@{
    //! Mesh typedefs
    typedef int                           dim_type;
    typedef size_t                        size_type;
    typedef profugus::geometry::cell_type cell_type;
    typedef cuda_utils::Space_Vector      Space_Vector;
    typedef cuda_utils::Coordinates       Coordinates;
    typedef std::vector<double>           Vec_Dbl;
    //@}

  private:
    // >>> DATA

    // Cell edges.
    thrust::device_vector<double> d_x_edges;
    thrust::device_vector<double> d_y_edges;
    thrust::device_vector<double> d_z_edges;

    // Dimensionality (always 3 for now)
    dim_type d_dimension;

    // Cell volumes
    std::vector<double> d_volumes;

    // Extents
    def::Space_Vector d_lower;
    def::Space_Vector d_upper;

  public:

    // Construct from xyz edges.
    Cartesian_Mesh_DMM(const Vec_Dbl& x_edges, const Vec_Dbl& y_edges,
                       const Vec_Dbl& z_edges);

    Cartesian_Mesh device_instance()
    {
        Cartesian_Mesh mesh(cuda::make_view(d_x_edges),
                            cuda::make_view(d_y_edges),
                            cuda::make_view(d_z_edges));
        return mesh;
    }

    // >>> ACCESSORS

    //! Number of cells in mesh
    size_type num_cells() const {return d_volumes.size();}

    // >>> VOLUME CALCULATION

    //! Get all volumes on host
    const std::vector<double>& volumes() const
    {
        return d_volumes;
    }

    //! Get lower corner of domain
    def::Space_Vector lower() const
    {
        return d_lower;
    }

    //! Get high corner of domain
    def::Space_Vector upper() const
    {
        return d_upper;
    }
};

} // end namespace cuda_profugus


#endif // cuda_geometry_Cartesian_Mesh_hh

//---------------------------------------------------------------------------//
//                 end of Cartesian_Mesh.hh
//---------------------------------------------------------------------------//
