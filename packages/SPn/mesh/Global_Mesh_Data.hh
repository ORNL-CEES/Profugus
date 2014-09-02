//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Global_Mesh_Data.hh
 * \author Thomas M. Evans
 * \date   Tue Feb 12 09:00:23 2014
 * \brief  Global_Mesh_Data class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mesh_Global_Mesh_Data_hh
#define mesh_Global_Mesh_Data_hh

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "Mesh.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Global_Mesh_Data
 * \brief Structure to store global mesh data.
 */
//===========================================================================//

class Global_Mesh_Data
{
  private:
    // >>> DATA

    // Global mesh edges.
    Mesh::Cell_Edges d_edges;

    // Low-edge boundaries in each dimension.
    Mesh::Space_Vector d_low;

    // High-edge boundaries in each dimension.
    Mesh::Space_Vector d_high;

    // Number of cells in each dimension.
    Mesh::Dim_Vector d_N;

    // Number of vertices in each dimension.
    Mesh::Dim_Vector d_V;

    // Total volume.
    double d_volume ;

    // dimensionality of problem (2 or 3)
    const int d_dimension;

  public:
    // >>> CONSTRUCTORS

    // 3D constructor
    inline Global_Mesh_Data(const def::Vec_Dbl &x_edges,
                            const def::Vec_Dbl &y_edges,
                            const def::Vec_Dbl &z_edges);

    // 2D constructor
    inline Global_Mesh_Data(const def::Vec_Dbl &x_edges,
                            const def::Vec_Dbl &y_edges);

    // >>> ACCESSORS

    //! Get dimension of problem ( 2 or 3 )
    int which_dimension() const { return d_dimension; }

    //! Get global mesh volume (cm\f$^{3}\f$).
    double volume() const { return d_volume; }

    //! Get the global mesh edges.
    const def::Vec_Dbl& edges(int dir) const { return d_edges[dir]; }

    //! Get global low-edge boundaries.
    double low_edge(int dir) const { return d_low[dir]; }

    //! Get global high-edge boundaries.
    double high_edge(int dir) const { return d_high[dir]; }

    //! Get global number of cells.
    int num_cells() const
    {
        REQUIRE(d_dimension == 2 || d_dimension == 3);
        int temp =1;
        for (int i=0; i<d_dimension;++i)
            temp *= d_N[i];
        return temp;
    }

    //! Get global number of vertices.
    int num_vertices() const
    {
        REQUIRE(d_dimension == 2 || d_dimension == 3);
        int temp=1;
        for (int i=0; i<d_dimension; ++i)
          temp *= d_V[i];
        return temp;
    }

    //! Get global number of cells in a direction.
    int num_cells(int dir) const { return d_N[dir]; }

    //! Get global number of vertices in a direction.
    int num_vertices(int dir) const { return d_V[dir]; }

}; //end Global_Mesh_data

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Global_Mesh_Data.i.hh"

#endif // mesh_Global_Mesh_Data_hh

//---------------------------------------------------------------------------//
//              end of mesh/Global_Mesh_Data.hh
//---------------------------------------------------------------------------//
