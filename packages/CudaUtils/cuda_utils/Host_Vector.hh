//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Host_Vector.hh
 * \author Seth R Johnson
 * \date   Mon Aug 12 08:48:53 2013
 * \brief  Host_Vector class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Host_Vector_hh
#define CudaUtils_cuda_utils_Host_Vector_hh

#include "harness/DBC.hh"
#include "Definitions.hh"

// Declare fields to avoid header propagation into CUDA kernel files
namespace profugus
{
template<typename T> class View_Field;
template<typename T> class const_View_Field;
}

namespace cuda
{
//===========================================================================//
namespace alloc
{
//! How to allocate host memory
enum Flag
{
    DEFAULT,         //!< Standard pinned memory
    MAPPED,          //!< No host/device copying is needed
    WRITE_COMBINED,  //!< Faster if you're just doing host->device
    MAPPED_WRITE_COMBINED //!< Mapped and write-combined
};
}

template<class Arch_T>
alloc::Flag adjust_alloc_flag(alloc::Flag);

//! For host emulation, don't use write-combined or mapped memory
template<>
inline alloc::Flag adjust_alloc_flag<cuda::arch::Host>(alloc::Flag)
{
    return alloc::DEFAULT;
}

//! For actual device mode, pass through the given flags
template<>
inline alloc::Flag adjust_alloc_flag<cuda::arch::Device>(alloc::Flag flag)
{
    return flag;
}

//===========================================================================//
/*!
 * \class Host_Vector
 * \brief Page-locked host memory container for faster CUDA transfers
 *
 * The benefits and drawbacks of the different allocation flags are described
 * in the CUDA manual. In a nutshell:
 *   - "mapped" memory on some platforms can obviate the copying of data to and
 *     from the GPU;
 *   - "write-combined" memory is faster to copy to the GPU but is write-only.
 *
 * On the host, you can access this data with iterators as well as the
 * assignment operator. On the device, if it uses mapped memory, you can acess
 * the GPU-local data with the same data() and cdata() accessors that
 * Device_Vector has.
 */
//===========================================================================//

template<class T>
class Host_Vector
{
    typedef Host_Vector<T> This;
  public:
    //@{
    //! Container typedefs
    typedef T        value_type;
    typedef T&       reference;
    typedef const T& const_reference;
    typedef T*       pointer;
    typedef const T* const_pointer;
    typedef T*       iterator;
    typedef const T* const_iterator;
    typedef size_t   size_type;
    //@}

    typedef profugus::const_View_Field<T> const_View_Field_t;
    typedef profugus::View_Field<T>       View_Field_t;
    typedef alloc::Flag                  Alloc_Flag;

  public:
    // Construct with number of elements (no initialization!)
    explicit Host_Vector(
            size_t count,
            const_reference value = T(),
            Alloc_Flag flag = alloc::DEFAULT);

    // Free memory on destruct
    ~Host_Vector();

    // Assign same-sized host memory
    void assign(const_View_Field_t hostvec);

    // Assign same-sized host memory from another Host vector
    void assign(const This& hostvec);

    // >>> ACCESSORS for host memory

    //! Field value at index i.
    reference operator[](size_type i)
    {
        REQUIRE(d_begin + i < d_end);
        return *(d_begin + i);
    }

    //! Constant field value at index i.
    const_reference operator[](size_type i) const
    {
        REQUIRE(!is_write_combined());
        REQUIRE(d_begin + i < d_end);
        return *(d_begin + i);
    }

    //! Number of field values.
    size_type size() const
    {
        return d_end - d_begin;
    }

    //! Whether the size of the field is zero
    bool empty() const { return d_begin == d_end; }

    //! Whether the data is available without copying to the GPU
    bool is_mapped() const
    {
        return (   d_alloc_flag == alloc::MAPPED
                || d_alloc_flag == alloc::MAPPED_WRITE_COMBINED);
    }

    //! Whether the data is available without copying to the GPU
    bool is_write_combined() const
    {
        return (   d_alloc_flag == alloc::WRITE_COMBINED
                || d_alloc_flag == alloc::MAPPED_WRITE_COMBINED);
    }

    //! Accessor to host memory location, for copying to CUDA only
    const_pointer cpu_data() const { return d_begin; }

  public:
    // >>> ACCESSORS for device if allocated as mapped memory

    // Get a device pointer if this is mapped memory
    pointer data()
    {
        REQUIRE(is_mapped());
        ENSURE(d_gpu_data);
        return d_gpu_data;
    }

    // Get a const device pointer if this is mapped memory
    const_pointer data() const
    {
        REQUIRE(is_mapped());
        ENSURE(d_gpu_data);
        return d_gpu_data;
    }

    // Get a const device pointer if this is mapped memory
    const_pointer cdata() const { return data(); }

  public:
    // >>> ITERATORS

    //! Begin iterator.
    iterator begin() { return d_begin; }

    //! Begin const_iterator.
    const_iterator begin() const
    {
        REQUIRE(!is_write_combined());
        return d_begin;
    }

    //! End iterator.
    iterator end() { return d_end; }

    //! End const_iterator.
    const_iterator end() const { return d_end; }

  private:
    // >>> DATA

    // Pointer to CUDA-allocated host memory
    pointer d_begin;

    // Pointer to begin + size
    pointer d_end;

    // Device pointer, if this is mapped memory
    pointer d_gpu_data;

    // Allocation type (used for checking whether we can get device pointer)
    const Alloc_Flag d_alloc_flag;

  private:
    // Disallow copying for now
    Host_Vector(const This&);
};

} // end namespace cuda

#endif // CudaUtils_cuda_utils_Host_Vector_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Host_Vector.hh
//---------------------------------------------------------------------------//
