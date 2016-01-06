//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_xs/test/XS_Device_Tester.hh
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_xs_test_XS_Device_Tester_hh
#define cuda_xs_test_XS_Device_Tester_hh

#include "cuda_utils/Shared_Device_Ptr.hh"
#include "xs/XS.hh"
#include "../cuda_xs/XS_Device.hh"

//---------------------------------------------------------------------------//
class XS_Device_Tester
{
  public:

    XS_Device_Tester( const profugus::XS& xs );

    const typename profugus::XS::Vector& vector( int matid, int type ) const;

    const typename profugus::XS::Matrix& matrix( int matid, int pn ) const;

    int pn_order() const { return d_xs.get_host_ptr()->pn_order(); }
    int num_groups() const { return d_xs.get_host_ptr()->num_groups(); }
    int num_mat() const { return d_xs.get_host_ptr()->num_mat(); }

  private:

    cuda::Shared_Device_Ptr<cuda_profugus::XS_Device> d_xs;
};

//---------------------------------------------------------------------------//

#endif // cuda_xs_test_XS_Device_Tester_hh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/XS_Device_Tester.hh
//---------------------------------------------------------------------------//
