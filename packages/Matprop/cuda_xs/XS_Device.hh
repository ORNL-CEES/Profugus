//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/XS_Device.hh
 * \author Stuart Slattery
 * \brief  XS class definition on the device.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef xs_XS_Device_hh
#define xs_XS_Device_hh

#include "cuda_utils/CudaMacros.hh"
#include "cuda_utils/SerialDenseDeviceVector.hh"
#include "cuda_utils/SerialDenseDeviceMatrix.hh"

#include "xs/XS.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class XS_Device
 * \brief Device-side Cross-section container class.
 */
/*!
 * \example xs/test/tstXS_Device.cc
 *
 * Test of XS_Device.
 */
//===========================================================================//
class XS_Device
{
  public:

    // Constructor.
    XS_Device( const profugus::XS& xs );

    // Destructor. Prohibits copy construction and assignment.
    ~XS_Device();

    // >>> ACCESSORS

    //! Pn order of data.
    PROFUGUS_HOST_DEVICE_FUNCTION
    int pn_order() const { return d_pn; }

    //! Number of groups in data.
    PROFUGUS_HOST_DEVICE_FUNCTION
    int num_groups() const { return d_Ng; }

    //! Number of materials in database.
    PROFUGUS_HOST_DEVICE_FUNCTION
    int num_mat() const { return d_Nm; }

    // Return the 1-D data vector for a given matid and type.
    PROFUGUS_DEVICE_FUNCTION
    const cuda::SerialDenseDeviceVector vector(int matid, int type) const
    { 
	PROFUGUS_INSIST_ON_DEVICE;
	int offset = type * d_Nm * d_Ng + d_matid_g2l[matid] * d_Ng;
	return cuda::SerialDenseDeviceVector( d_Ng, d_totals+offset );
    }

    // Return the 2-D data matrix for a given matid and Pn order.
    PROFUGUS_DEVICE_FUNCTION
    const cuda::SerialDenseDeviceMatrix matrix(int matid, int pn) const
    { 
	PROFUGUS_INSIST_ON_DEVICE;
	int offset = pn * d_Nm * d_Ng * d_Ng + d_matid_g2l[matid] * d_Ng * d_Ng;
	return cuda::SerialDenseDeviceMatrix( d_Ng, d_Ng, d_scatter+offset );
    }

  private:

    // >>> DATA

    // Pn order of scattering.
    int d_pn;

    // Number of groups.
    int d_Ng;

    // Final number of materials.
    int d_Nm;

    // Number of cross section types.
    int d_Nxst;

    // Global-to-local mapping of matids. On-device.
    int* d_matid_g2l;

    // Size of total cross section data.
    std::size_t d_totals_size;

    // Total cross sections on device strided by type and group size. On-device.
    double* d_totals;

    // Size of scattering cross section data.
    std::size_t d_scatter_size;

    // Scattering cross sections on device. On-device.
    double* d_scatter;

};

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // xs_XS_Device_hh

//---------------------------------------------------------------------------//
//                 end of XS_Device.hh
//---------------------------------------------------------------------------//
