//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Mesh.hh
 * \author Tom Evans
 * \date   Wednesday February 12 0:15:45 2014
 * \brief  Declaration of Mesh
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Mesh_kba_hh
#define spn_Mesh_kba_hh

#include <vector>

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "utils/Vector_Lite.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Mesh
 * \brief Orthogonal mesh class.
 *
 * The mesh is partitioned in (I,J) blocks where I has range [0,Np_I) and J
 * has range [0,Np_J) and the \f$ N_p = (N_p)_I\times(N_p)_J\f$.  The
 * processor id for this mesh is calculated using
 * \f[
 * \mbox{proc\_id} =
 * \mbox{I\_block} + \mbox{J\_block}\times(N_p)_I\:.
 * \f]
 *
 * There are no separate domains in K; however, the user can define effective
 * blocks in K for added parallel efficiency.  The number of cells in each
 * "effective" parallel block is thus \f$ N_x\times N_y\times\mbox{int}(N_z /
 * \mbox{num\_K\_blocks})\f$.  All meshes in the decomposition must have the
 * same number of K blocks and \f$ N_z \f$ defined.  A current restriction is
 * that \f$ N_z\, \mbox{mod}\,\mbox{num\_K\_blocks} = 0 \f$.
 *
 */
/*!
 * \example spn/test/tstMesh.cc
 *
 * Test of Mesh.
 */
//===========================================================================//

class Mesh
{
  public:

    //@{
    //! Useful typedefs for classes that use mesh.
    typedef def::size_type size_type;
    typedef def::tiny_int  dim_type;

    typedef profugus::Vector_Lite<int, 3>          Dim_Vector;
    typedef profugus::Vector_Lite<double, 3>       Space_Vector;
    typedef profugus::Vector_Lite<def::Vec_Dbl, 3> Cell_Edges;
    typedef profugus::Vector_Lite<def::Vec_Dbl, 3> Cell_Widths;
    //@}

  public:
    // Build a 3D mesh (or, if z_edges is empty, a 2-D mesh).
    Mesh(const def::Vec_Dbl &x_edges,
         const def::Vec_Dbl &y_edges,
         const def::Vec_Dbl &z_edges,
         size_type           I_block,
         size_type           J_block,
         size_type           num_K_blocks);

    // Build a 2D mesh
    Mesh(const def::Vec_Dbl &x_edges,
         const def::Vec_Dbl &y_edges,
         size_type           I_block,
         size_type           J_block);

    // >>> PUBLIC INTERFACE

    // Convert cardinal index to (i,j) or (i,j,k).
    inline Dim_Vector cardinal(size_type cell) const;

    // Cell volume.
    inline double volume(size_type cell) const;

    // Block width.
    inline Space_Vector block_widths() const;

    // Query to see if a point is in a mesh.
    bool find_cell(const Space_Vector &r, Dim_Vector &ijk) const;

    // Label of mesh
    std::string label() const
    {
        if (dimension() == 3)
            return "xyz_3d";
        return "xy_2d";
    }

    // >>> NON-VIRTUAL FUNCTIONS

    //! Dimension of mesh.
    size_type dimension() const { return d_dimension; }

    //! Get number of cells.
    size_type num_cells() const { return d_num_cells; }

    //! Get number of vertices.
    size_type num_vertices() const
    {
        REQUIRE(d_num_vertices != 0);
        return d_num_vertices;
    }

    // Returns the width of a cell in a given direction
    inline double width(size_type ijk, int dim) const;

    // Return inverse of cell width along \e (i,j) or \e (i,j,k).
    inline double inv_width(size_type ijk, int dim) const;

    // convert (i,j,k) to cardinal index
    inline size_type convert(size_type i, size_type j, size_type k) const;

    // Convert (i,j) to cardinal index assumes k=0 (if 3D)
    inline size_type convert(size_type i, size_type j) const;

    // Cell center along \e (i,j) or \e(i,j,k).
    inline double center(size_type ijk, int dim) const;

    //! Low corner of mesh in \e (i,j,k) direction.
    double low_corner(int d) const
    {
        REQUIRE(d < static_cast<int>(d_dimension) ); return d_edges[d].front();
    }

    //! High corner of mesh in \e (i,j,k) direction.
    double high_corner(int d) const
    {
        REQUIRE(d < static_cast<int>(d_dimension) ); return d_edges[d].back();
    }

    //! Return all cell edges.
    const Cell_Edges& edges() const { return d_edges; }

    //! Return cell edges along a given direction.
    const def::Vec_Dbl& edges(int d) const
    {
        REQUIRE(d < static_cast<int>(d_dimension) ); return d_edges[d];
    }

    //! Return cell edge coordinate along \e d.
    double edges(int vertices, int d) const
    {
        REQUIRE(d < static_cast<int>(d_dimension) &&
                vertices < static_cast<int>(d_edges[d].size()));
        return d_edges[d][vertices];
    }

    //! Return block width in a particular direction.
    double block_width(int d) const
    {
        REQUIRE(d < 3);
        return d_edges[d].back() - d_edges[d].front();
    }

    //! Return number of cells in each dimension in this mesh.
    const Dim_Vector& num_cells_dims() const { return d_N; }

    //! Return number of cells in a given dimension \e ijk in this mesh.
    inline int num_cells_dim(size_type ijk) const
    {
        REQUIRE(ijk < 3); return d_N[ijk];
    }

    //! Get block description (block index on i/j; num blocks along z).
    const Dim_Vector& get_blocks() const { return d_blocks; }

    //! Get block index along dimension (or number of blocks if 3)
    size_type block(size_type ijk) const
    {
        REQUIRE(ijk < d_dimension); return d_blocks[ijk];
    }

    //! Get number of cells in each dimension in each block in this mesh.
    const Dim_Vector& num_cells_block_dims() const { return d_Nb; }

    //! Get number of cells in a given block along dimension \e (i,j,k).
    size_type num_cells_block_dim(size_type ijk) const
    {
        REQUIRE(ijk < d_dimension); return d_Nb[ijk];
    }

  protected:
    // >>> protected interfaces used in the derived mesh classes

    // builds cell widths and inverse cell widths
    void build_mesh();

  private:
    // >>> DATA

    // dimension of mesh
    const size_type d_dimension;

    // Number of cells in mesh.
    size_type d_num_cells;

    // Number of vertices in mesh.
    size_type d_num_vertices;

    //! Number of cells in each dimension in this mesh.
    Dim_Vector d_N;

    //@{
    //! Block description either IxJ if 2D or IxJxK in 3D, I J are block IDs,
    //! K is the number of blocks in the Z direction
    Dim_Vector d_blocks;
    //@}

    //! Number of cells in each dimension in each block of this mesh.
    Dim_Vector d_Nb;

    // Cell edges.
    Cell_Edges d_edges;

    // Cell widths.
    Cell_Widths d_width;

    // Inverse cell widths.
    Cell_Widths d_inv_width;
};

//---------------------------------------------------------------------------//

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Mesh.i.hh"

#endif // spn_Mesh_hh

//---------------------------------------------------------------------------//
//              end of spn/Mesh.hh
//---------------------------------------------------------------------------//
