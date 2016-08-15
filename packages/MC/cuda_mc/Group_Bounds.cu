//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Group_Bounds.cu
 * \author Stuart Slattery
 * \brief  Group_Bounds member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <algorithm>

#include "Group_Bounds.hh"

#include "cuda_utils/Memory.cuh"

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
/*!
 * \brief Construct with neutron and photon boundaries.
 *
 * Boundaries must be monotonic decreasing or empty. At least one group
 * boundaries must be provided.
 */
Group_Bounds::Group_Bounds(const Vec_Dbl &bounds)
    : d_size( bounds.size() )
{
    // check monotonicity
    for (int g = 0; g < num_groups(); ++g)
    {
        CHECK( bounds[g] > bounds[g+1] );
    }

    // Copy to device.
    cuda_utils::memory::Malloc( d_bounds, bounds.size() );
    cuda_utils::memory::Copy_To_Device( d_bounds, bounds.data(), bounds.size() );
}

//---------------------------------------------------------------------------//
/*
 * \brief Destructor.
 */
Group_Bounds::~Group_Bounds()
{
    cuda_utils::memory::Free( d_bounds );
}

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
//                 end of Group_Bounds.cu
//---------------------------------------------------------------------------//
