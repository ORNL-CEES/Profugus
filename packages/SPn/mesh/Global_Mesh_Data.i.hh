//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/mesh/Global_Mesh_Data.i.hh
 * \author Thomas M. Evans
 * \date   Tue Feb 12 09:00:23 2008
 * \brief  Global_Mesh_Data class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_mesh_Global_Mesh_Data_i_hh
#define SPn_mesh_Global_Mesh_Data_i_hh

namespace profugus
{

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief 3D Constructor
 */
Global_Mesh_Data::Global_Mesh_Data(const def::Vec_Dbl &x_edges,
                                   const def::Vec_Dbl &y_edges,
                                   const def::Vec_Dbl &z_edges)
    : d_edges  ( x_edges, y_edges, z_edges)
    , d_low    ( x_edges.front(), y_edges.front(), z_edges.front() )
    , d_high   ( x_edges.back() , y_edges.back() , z_edges.back()  )
    , d_N      ( x_edges.size() - 1 , y_edges.size() - 1 , z_edges.size() - 1)
    , d_V      ( d_N[0] + 1 , d_N[1] + 1 , d_N[2] + 1 )
    , d_volume ( (d_high[0]-d_low[0]) *
                 (d_high[1]-d_low[1]) *
                 (d_high[2]-d_low[2]) )
    , d_dimension(3)
{
} //end constructor 3D

//---------------------------------------------------------------------------//
/*!
 * \brief 2D Constructor
 */
Global_Mesh_Data::Global_Mesh_Data(const def::Vec_Dbl &x_edges,
                                   const def::Vec_Dbl &y_edges)
    : d_edges  ( x_edges, y_edges, def::Vec_Dbl(0) )
    , d_low    ( x_edges.front(), y_edges.front(), 0.0 )
    , d_high   ( x_edges.back() , y_edges.back() , 0.0 )
    , d_N      ( x_edges.size() - 1 , y_edges.size() - 1 , 0 )
    , d_V      ( d_N[0] + 1 , d_N[1] + 1 , 0 )
    , d_volume ( (d_high[0]-d_low[0]) * (d_high[1]-d_low[1]) )
    , d_dimension(2)
{
} //end constructor 3D

} // end namespace profugus

#endif // SPn_mesh_Global_Mesh_Data_i_hh

//---------------------------------------------------------------------------//
//              end of mesh/Global_Mesh_Data.i.hh
//---------------------------------------------------------------------------//
