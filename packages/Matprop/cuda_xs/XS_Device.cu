//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/XS_Device.cu
 * \author Stuart Slattery
 * \brief  XS class definition on the device.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "XS_Device.hh"

#include "cuda_utils/Memory.cuh"

#include <Teuchos_Array.hpp>

#include <cuda_runtime.h>

#include <algorithm>

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
XS_Device::XS_Device( const profugus::XS& xs )
    : d_pn( xs.pn_order() )
    , d_Ng( xs.num_groups() )
    , d_Nm( xs.num_mat() )
    , d_Nxst( profugus::XS::END_XS_TYPES )
    , d_totals_size( d_Nxst * d_Nm * d_Ng )
    , d_scatter_size( (d_pn+1) * d_Nm * d_Ng * d_Ng)
{
    // Get the matids.
    typename profugus::XS::Vec_Int matids;
    xs.get_matids( matids );

    // Create a global to local mapping of matids.
    int matid_g2l_size = *std::max_element( matids.begin(), matids.end() ) + 1;
    Teuchos::Array<int> host_matid_g2l( matid_g2l_size, -1 );
    for ( int m = 0; m < d_Nm; ++m )
    {
	host_matid_g2l[ matids[m] ] = m;
    }

    // Allocate a matid global-to-local map.
    cuda::memory::Malloc( d_matid_g2l, matid_g2l_size );

    // Copy the matid list to the device.
    cuda::memory::Copy_To_Device( 
	d_matid_g2l, host_matid_g2l.getRawPtr(), matid_g2l_size );
    host_matid_g2l.clear();

    // Allocate total cross sections.
    cuda::memory::Malloc( d_totals, d_totals_size );

    // Extract the total cross sections.
    double* host_xs;
    std::size_t offset = 0;
    for ( int t = 0; t < d_Nxst; ++t )
    {
	for ( int m = 0; m < d_Nm; ++m )
	{
	    host_xs = xs.vector( matids[m], t ).values();
	    offset = t * d_Nm * d_Ng + m * d_Ng;
	    cuda::memory::Copy_To_Device( d_totals + offset, host_xs, d_Ng );
	}
    }

    // Allocate the scattering cross sections.
    cuda::memory::Malloc( d_scatter, d_scatter_size );

    // Extract the scattering cross sections.
    for ( int pn = 0; pn < d_pn+1; ++pn )
    {
	for ( int m = 0; m < d_Nm; ++m )
	{
	    host_xs = xs.matrix( matids[m], pn ).values();
	    offset = pn * d_Nm * d_Ng * d_Ng + m * d_Ng * d_Ng;
	    cuda::memory::Copy_To_Device( 
		d_scatter + offset, host_xs, d_Ng*d_Ng );
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
XS_Device::~XS_Device()
{
    cuda::memory::Free( d_matid_g2l );
    cuda::memory::Free( d_totals );
    cuda::memory::Free( d_scatter );
}

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
//                 end of XS_Device.hh
//---------------------------------------------------------------------------//
