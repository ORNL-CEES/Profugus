//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/Manager_Builder_Cuda.cu
 * \author Steven Hamilton
 * \date   Wed Mar 30 15:02:55 2016
 * \brief  Manager_Builder_Cuda class definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Manager_Builder_Cuda.hh"
#include "Manager_Cuda.hh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_mc
{
//---------------------------------------------------------------------------//

auto Manager_Builder_Cuda::build( const RCP_ParameterList &pl )
    -> SP_Manager_Base
{
    typedef cuda_profugus::Mesh_Geometry Mesh;

    if( pl->isSublist("MESH") )
        return std::make_shared<Manager_Cuda<Mesh>>();
    else
        INSIST(false,"Only MESH geometry available in cuda.");

    return SP_Manager_Base();
}

//---------------------------------------------------------------------------//
} // end namespace cuda_mc

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Manager_Builder_Cuda.cu
//---------------------------------------------------------------------------//
