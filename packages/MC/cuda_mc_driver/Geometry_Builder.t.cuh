//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc_driver/Geometry_Builder.t.cuh
 * \author Steven Hamilton
 * \date   Wed Nov 25 12:58:58 2015
 * \brief  Geometry_Builder template member definitions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_driver_Geometry_Builder_t_hh
#define cuda_mc_driver_Geometry_Builder_t_hh

#include "Geometry_Builder.hh"

namespace cuda_mc
{

auto Geometry_Builder<cuda_profugus::Mesh_Geometry>::build(
    RCP_ParameterList master) -> SDP_Geometry
{
    auto mesh_db = Teuchos::sublist(master, "MESH");

    // Ensure all required parameters are present
    REQUIRE( mesh_db->isParameter("x_edges") );
    REQUIRE( mesh_db->isParameter("y_edges") );
    REQUIRE( mesh_db->isParameter("z_edges") );
    REQUIRE( mesh_db->isParameter("matids") );

    auto x_edges = mesh_db->get<OneDArray_dbl>("x_edges");
    auto y_edges = mesh_db->get<OneDArray_dbl>("y_edges");
    auto z_edges = mesh_db->get<OneDArray_dbl>("z_edges");
    auto matids  = mesh_db->get<OneDArray_int>("matids");

    // Build Mesh
    auto geom = std::make_shared<cuda_profugus::Mesh_Geometry>(
        x_edges.toVector(),y_edges.toVector(),z_edges.toVector());

    // Convert matids to SP<Vec_Int> and pass to geometry
    auto sp_matids = std::make_shared<def::Vec_Int>(
        matids.begin(),matids.end());

    REQUIRE( sp_matids->size() == ( (x_edges.size()-1) *
                                    (y_edges.size()-1) *
                                    (z_edges.size()-1) ) );
    geom->set_matids(sp_matids);

    return SDP_Geometry(geom);
}

} // end namespace cuda_mc

#endif // cuda_mc_driver_Geometry_Builder_t_hh

//---------------------------------------------------------------------------//
//                 end of Geometry_Builder.t.cuh
//---------------------------------------------------------------------------//
