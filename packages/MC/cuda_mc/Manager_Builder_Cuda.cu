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
#include "cuda_rtk/RTK_Geometry.cuh"

namespace cuda_mc
{
//---------------------------------------------------------------------------//

auto Manager_Builder_Cuda::build( const RCP_ParameterList &pl )
    -> SP_Manager_Base
{
    typedef cuda_profugus::Mesh_Geometry_DMM Mesh_DMM;
    typedef cuda_profugus::Core_DMM          Core_DMM;

    if (pl->isSublist("MESH"))
        return std::make_shared<Manager_Cuda<Mesh_DMM>>();
    else if (pl->isSublist("CORE"))
        return std::make_shared<Manager_Cuda<Core_DMM>>();
    else
        VALIDATE(false,"Unrecognized geometry type.");

    return SP_Manager_Base();
}

//---------------------------------------------------------------------------//
} // end namespace cuda_mc

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Manager_Builder_Cuda.cu
//---------------------------------------------------------------------------//
