//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/Cartesian_Mesh.hh
 * \author Thomas M. Evans
 * \date   Mon Jul 21 16:55:00 2014
 * \brief  Cartesian_Mesh class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_Cartesian_Mesh_hh
#define geometry_Cartesian_Mesh_hh

#include <vector>
#include <cstddef>

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "utils/Vector_Lite.hh"

namespace profugus
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
    typedef profugus::Vector_Lite<dim_type, 3> Dim_Vector;
    typedef profugus::Vector_Lite<double, 3>   Space_Vector;
    typedef profugus::Vector_Lite<Vec_Dbl, 3>  Cell_Edges;
    //@}

  private:
    // >>> DATA

    // Cell edges.
    Cell_Edges d_edges;

    // Number of cells along each axis
    Dim_Vector d_extents;

    // Number of global cells
    size_type d_num_cells;

    // Dimensionality
    dim_type d_dimension;

  public:
    // Construct from xyz edges.
    Cartesian_Mesh(const Vec_Dbl& x_edges, const Vec_Dbl& y_edges,
                   const Vec_Dbl& z_edges);

    // >>> ACCESSORS

    //! Get number of cells.
    size_type num_cells() const { return d_num_cells; }

    //! Number of cells along an axis
    dim_type num_cells_along(dim_type d) const
    {
        REQUIRE(0 <= d && d < dimension());
        return d_extents[d];
    }

    //! Return number of cells in each dimension in this mesh.
    const Dim_Vector& extents() const { return d_extents; }

    //! Dimension of mesh.
    dim_type dimension() const { return d_dimension; }

    //! Return all cell edges.
    const Cell_Edges& edges() const { return d_edges; }

    //! Return cell edges along a given direction.
    const Vec_Dbl& edges(dim_type d) const
    {
        REQUIRE(0 <= d && d < d_dimension);
        return d_edges[d];
    }

    // >>> INDEX CONVERSION

    // Convert cardinal index to (i,j) or (i,j,k).
    void cardinal(size_type cell, dim_type& i, dim_type& j, dim_type& k) const;

    // Convert cardinal index to (i,j) or (i,j,k).
    inline Dim_Vector cardinal(size_type cell) const;

    // Convert (i,j,k) to cell index
    inline bool index(dim_type i, dim_type j, dim_type k,
                      size_type& cell) const;

    // Convert (i,j,k) to cell index, raise on bad index
    inline size_type index(dim_type i, dim_type j, dim_type k) const;

    // >>> VOLUME CALCULATION

    // Calculate cell volume.
    inline double volume(dim_type i, dim_type j, dim_type k) const;

    //! Calculate volume from the global cell id
    inline double volume(size_type global_cell_index) const;

    // >>> SPATIAL LOCATION

    // Locate a coordinate along a single axis
    inline dim_type find_upper(double r, dim_type axis) const;

    // Locate the position's ijk coordinates with upper edges being "inside"
    void find_upper(const Space_Vector& r, Dim_Vector& ijk) const;

    // Locate the position's ijk coordinates
    bool find(const Space_Vector& r, Dim_Vector& ijk) const;

    //! Low corner of mesh in \e (i,j,k) direction.
    double low_corner(dim_type d) const
    {
        REQUIRE(0 <= d && d < d_dimension);
        return d_edges[d].front();
    }

    //! High corner of mesh in \e (i,j,k) direction.
    double high_corner(dim_type d) const
    {
        REQUIRE(0 <= d && d < d_dimension);
        return d_edges[d].back();
    }
};

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Cartesian_Mesh.i.hh"

#endif // geometry_Cartesian_Mesh_hh

//---------------------------------------------------------------------------//
//                 end of Cartesian_Mesh.hh
//---------------------------------------------------------------------------//
