//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Device_Memory_Manager.hh
 * \author Tom Evans
 * \date   Thu Nov 17 15:50:00 2016
 * \brief  Device_Memory_Manager class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Device_Memory_Manager_hh
#define CudaUtils_cuda_utils_Device_Memory_Manager_hh

#include <type_traits>

namespace cuda_utils
{

//===========================================================================//
/*!
 * \class Device_Handle
 * \brief Opaque host handle to device memory.
 *
 * This is a base class that provides an opaque handle, stored on the host, to
 * device memory.  It also manages device object construction.  This
 * guarantees that deep objects are correctly deleted.
 */
//===========================================================================//

template<class T>
class Device_Memory_Manager
{
  public:
    //! Constructor.
    Device_Memory_Manager()
    {
        // NVCC does not support is_trivially_copyable<T>, but is ok with this
        // macro, super
        static_assert(__has_trivial_copy(T),
                      "T is not trivially copyable");
    }

    //! Destructor.
    virtual ~Device_Memory_Manager() = default;

    //! Interface.
    virtual T device_instance() = 0;
};

//---------------------------------------------------------------------------//

} // end namespace cuda_utils

#endif // CudaUtils_cuda_utils_Device_Memory_Manager_hh

//---------------------------------------------------------------------------//
// end of CudaUtils/cuda_utils/Device_Memory_Manager.hh
//---------------------------------------------------------------------------//
