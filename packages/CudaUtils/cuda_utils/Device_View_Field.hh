//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Device_View_Field.hh
 * \author Steven Hamilton, Tom Evans
 * \date   Thu Nov 17 16:30:59 2016
 * \brief  Device_View_Field class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Device_View_Field_hh
#define CudaUtils_cuda_utils_Device_View_Field_hh

#include "thrust/device_vector.h"

#include "Utils/harness/DBC.hh"
#include "CudaDBC.hh"

namespace cuda_utils
{

//===========================================================================//
/*!
 * \class Device_View_Field
 * \brief Non-owning view into contiguous device-allocated data.
 *
 * Note that unlike nemesis::View_Field, the Device_View_Field does not
 * support strided data.  This class is intended as a lightweight object
 * to encapsulate device data.  The pointers must point to device-allocated
 * data (we cannot verify this).  Pointers/iterators may be accessed from
 * either host or device, but can only be dereferenced on device.
 */
/*!
 * \example cuda_utils/test/tstDevice_View_Field.cc
 *
 * Test of Device_View_Field.
 */
//===========================================================================//

template <typename T>
class Device_View_Field
{
  public:
    //@{
    //! Typedefs
    typedef T           value_type;
    typedef T&          reference;
    typedef const T&    const_reference;
    typedef T*          pointer;
    typedef const T*    const_pointer;
    typedef int         size_type;

    typedef Device_View_Field<T>    Device_View_Field_t;
    //@}

  private:
    // >>> DATA
    pointer d_begin;
    pointer d_end;

  public:

    // Default Constructor
    __host__ __device__ Device_View_Field()
      : d_begin(nullptr)
      , d_end(nullptr)
    {}

    // Constructor from pointer and size
    __host__ __device__ Device_View_Field(pointer p, size_type size)
      : d_begin(p)
      , d_end(p+size)
    {
        DEVICE_REQUIRE(size >= 0);
    }

    //>>> Container functions

    // Element access
    __device__ reference operator[](size_type i)
    {
        DEVICE_REQUIRE(valid_index(i));
        return d_begin[i];
    }

    // Element access
    __device__ const_reference operator[](size_type i) const
    {
        DEVICE_REQUIRE(valid_index(i));
        return d_begin[i];
    }

    // First element
    __device__ reference front()
    {
        DEVICE_REQUIRE(!empty());
        return *d_begin;
    }

    // First element
    __device__ const_reference front() const
    {
        DEVICE_REQUIRE(!empty());
        return *d_begin;
    }

    // Last element
    __device__ reference back()
    {
        DEVICE_REQUIRE(!empty());
        return *(d_end-1);
    }

    // First element
    __device__ const_reference back() const
    {
        DEVICE_REQUIRE(!empty());
        return *(d_end-1);
    }

    // Size
    __host__ __device__ size_type size() const
    {
        return d_end - d_begin;
    }

    // Is container empty?  Device-only because this is called from accessors.
    __device__ bool empty() const
    {
        return d_begin == d_end;
    }

    // Begin iterator
    __host__ __device__ pointer       begin()       {return d_begin;}
    __host__ __device__ const_pointer begin() const {return d_begin;}

    // End iterator
    __host__ __device__ pointer       end()       {return d_end;}
    __host__ __device__ const_pointer end() const {return d_end;}

  private:

    __device__ bool valid_index(size_type i) const
    {
        return d_begin + i < d_end;
    }

};

//===========================================================================//
/*!
 * \class const_Device_View_Field
 * \brief Non-owning constant view into contiguous device-allocated data.
 *
 * Note that unlike nemesis::const_View_Field, the const_Device_View_Field
 * does not support strided data.  This class is intended as a lightweight
 * object to encapsulate device data.  The pointers must point to
 * device-allocated data (we cannot verify this).  Pointers/iterators may be
 * accessed from either host or device, but can only be dereferenced on device.
 */
//===========================================================================//

template <typename T>
class const_Device_View_Field
{
  public:
    //@{
    //! Typedefs
    typedef T           value_type;
    typedef T&          reference;
    typedef const T&    const_reference;
    typedef T*          pointer;
    typedef const T*    const_pointer;
    typedef int         size_type;

    typedef Device_View_Field<T>       Device_View_Field_t;
    typedef const_Device_View_Field<T> const_Device_View_Field_t;
    //@}

  private:
    // >>> DATA
    const_pointer d_begin;
    const_pointer d_end;

  public:

    // Default Constructor
    __host__ __device__ const_Device_View_Field()
      : d_begin(nullptr)
      , d_end(nullptr)
    {}

    // Constructor from a Device_View_Field
    __host__ __device__ const_Device_View_Field(const Device_View_Field_t &view)
      : d_begin(view.begin())
      , d_end(view.end())
    {}

    // Constructor from pointer and size
    __host__ __device__ const_Device_View_Field(const_pointer p, size_type size)
      : d_begin(p)
      , d_end(p+size)
    {
        DEVICE_REQUIRE(size >= 0);
    }

    // Assignment from another const_Device_View_Field
    __host__ __device__ const_Device_View_Field_t & operator=(
            const const_Device_View_Field_t &view)
    {
        d_begin = view.d_begin;
        d_end   = view.d_end;
        return *this;
    }

    // Assignment from another const_Device_View_Field
    __host__ __device__ const_Device_View_Field_t & operator=(
            const Device_View_Field_t &view)
    {
        d_begin = view.begin();
        d_end   = view.end();
        return *this;
    }

    //>>> Container functions

    // Element access
    __device__ const_reference operator[](size_type i) const
    {
        DEVICE_REQUIRE(valid_index(i));
        return d_begin[i];
    }

    // First element
    __device__ const_reference front() const
    {
        DEVICE_REQUIRE(!empty());
        return *d_begin;
    }

    // First element
    __device__ const_reference back() const
    {
        DEVICE_REQUIRE(!empty());
        return *(d_end-1);
    }

    // Size
    __host__ __device__ size_type size() const
    {
        return d_end - d_begin;
    }

    // Is container empty?
    __host__ __device__ bool empty() const
    {
        return d_begin == d_end;
    }

    // Begin iterator
    __host__ __device__ const_pointer begin() const {return d_begin;}

    // End iterator
    __host__ __device__ const_pointer end() const {return d_end;}

  private:

    __device__ bool valid_index(size_type i) const
    {
        return d_begin + i < d_end;
    }

};

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// Make Device_View_Field from thrust::device_vector
template <typename T>
Device_View_Field<T> make_view(thrust::device_vector<T> &vec)
{
    return Device_View_Field<T>(vec.data().get(),vec.size());
}

// Make const_Device_View_Field from thrust::device_vector
template <typename T>
const_Device_View_Field<T> make_view(const thrust::device_vector<T> &vec)
{
    return const_Device_View_Field<T>(vec.data().get(),vec.size());
}

// Make const_Device_View_Field from thrust::device_vector
template <typename T>
const_Device_View_Field<T> make_const_view(const thrust::device_vector<T> &vec)
{
    return const_Device_View_Field<T>(vec.data().get(),vec.size());
}

//---------------------------------------------------------------------------//

} // end namespace cuda_utils

#endif // CudaUtils_cuda_utils_Device_View_Field_hh

//---------------------------------------------------------------------------//
// end of CudaUtils/cuda_utils/Device_View_Field.hh
//---------------------------------------------------------------------------//
