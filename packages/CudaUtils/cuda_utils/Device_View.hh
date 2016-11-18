//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Device_View.hh
 * \author Tom Evans
 * \date   Thu Nov 17 16:21:25 2016
 * \brief  Device_View class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Device_View_hh
#define CudaUtils_cuda_utils_Device_View_hh

#include <memory>
#include <vector>
#include <type_traits>

#include <thrust/device_vector.h>

#include "harness/DBC.hh"
#include "Device_Memory_Manager.hh"
#include "Device_View_Field.hh"

namespace cuda
{

//===========================================================================//
/*!
 * \class Device_View
 * \brief View to device objects.
 *
 * This class allows the client to construct a field of device objects.  It
 * uses thrust to manage memory on the device.
 */
/*!
 * \example cuda_utils/test/tstDevice_View.cc
 *
 * Test of Device_View.
 */
//===========================================================================//

template<class T>
class Device_View
{
  public:
    //@{
    //! Typedefs.
    using SP_Memory_Manager = std::shared_ptr<Device_Memory_Manager<T>>;
    using Managers          = std::vector<SP_Memory_Manager>;
    using Field             = std::vector<T>;
    using const_pointer     = typename Field::const_pointer;
    using pointer           = typename Field::pointer;
    using size_type         = typename Field::size_type;
    //@}

  private:
    // >>> DATA

    //! Managers.
    Managers d_managers;

    //! Field data.
    thrust::device_vector<T> d_field;

  public:
    // Constructor
    Device_View(Managers managers)
        : d_managers(std::move(managers))
    {
        Field field;
        for (auto &manager : d_managers)
        {
            field.push_back(manager->device_instance());
        }

        d_field = thrust::device_vector<T>(field.begin(), field.end());
    }

    // >>> ACCESSORS

    // Get Device_View_Field to data
    Device_View_Field<T>       get_view()
    {
        return cuda::make_view(d_field);
    }

    const_Device_View_Field<T> get_view() const
    {
        return cuda::make_view(d_field);
    }

    //@{
    //! Return pointer to device field.
    const_pointer get_device_ptr() const { return d_field.data().get(); }
    pointer get_device_ptr() { return d_field.data().get(); }
    //@}

    //! Get the size of the field.
    size_type size() const { return d_field.size(); }
};

//---------------------------------------------------------------------------//

} // end namespace cuda

#endif // CudaUtils_cuda_utils_Device_View_hh

//---------------------------------------------------------------------------//
// end of CudaUtils/cuda_utils/Device_View.hh
//---------------------------------------------------------------------------//
