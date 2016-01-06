//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/XS_Device.cu
 * \author Stuart Slattery
 * \brief  XS class definition on the device.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "XS_Device.hh"

#include <Teuchos_Array.hpp>

#include <cuda_runtime.h>

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
    , d_scatter_size( d_pn * d_Nm * d_Ng * d_Ng)
{
    // Get the matids.
    typename profugus::XS::Vec_Int matids;
    xs.get_matids( matids );

    // Create a global to local mapping of matids.
    Teuchos::Array<int> host_matid_g2l( d_Nm, -1 );
    for ( int m = 0; i < d_Nm; ++m )
    {
	host_matid_g2l[ matids[m] ] = m;
    }

    // Allocate a matid global-to-local map.
    cudaMalloc( (void**) &d_matid_g2l, d_Nm*sizeof(int) );

    // Copy the matid list to the device.
    cudaMemcpy( d_matid_g2l, host_matid_g2l, d_Nm*sizeof(int),
		cudaMemcpyHostToDevice );
    host_matid_g2l.clear();

    // Allocate total cross sections.
    cudaMalloc( (void**) &d_totals, d_totals_size*sizeof(double) );

    // Extract the total cross sections.
    double* host_xs;
    std::size_t offset = 0;
    for ( int t = 0; t < d_Nxst; ++t )
    {
	for ( int m = 0; m < d_Nm; ++m )
	{
	    host_xs = xs.vector( matids[m], t ).values();
	    offset = t * d_Nm * d_Ng + m * d_Ng;
	    cudaMemcpy( d_totals + offset, host_xs, d_Ng*sizeof(double),
			cudaMemcpyHostToDevice );
	}
    }

    // Allocate the scattering cross sections.
    cudaMalloc( (void**) &d_scatter, d_scatter_size*sizeof(double) );

    // Extract the scattering cross sections.
    for ( int pn = 0; pn < d_pn; ++pn )
    {
	for ( int m = 0; m < d_Nm; ++m )
	{
	    host_xs = xs.matrix( matids[m], pn ).values();
	    offset = pn * d_Nm * d_Ng * d_Ng + m * d_Ng * d_Ng;
	    cudaMemcpy( d_scatter + offset, host_xs, d_Ng*d_Ng*sizeof(double),
			cudaMemcpyHostToDevice );
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
XS_Device::~XS_Device()
{
    cudaFree( d_matid_g2l );
    cudaFree( d_totals );
    cudaFree( d_scatter );
}

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
//                 end of XS_Device.hh
//---------------------------------------------------------------------------//
