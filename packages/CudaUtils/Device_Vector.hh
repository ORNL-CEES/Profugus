//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Device_Vector.hh
 * \author Seth R Johnson
 * \date   Thu Aug 01 11:33:12 2013
 * \brief  Device_Vector class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Device_Vector_hh
#define cuda_utils_Device_Vector_hh

#include <cuda_utils/config.h>

#ifndef __NVCC__
#include <vector>
#endif

#include "Utils/harness/DBC.hh"

#include "Definitions.hh"

// Declare fields to avoid header propagation into CUDA kernel files
namespace profugus
{
template<typename T> class View_Field;
template<typename T> class const_View_Field;
}

namespace cuda
{
// Declare host vector
template <typename T> class Host_Vector;
// Declare stream
template <typename Arch_T> class Stream;

//===========================================================================//
/*!
 * \class Device_Vector
 * \brief Heap-allocated vector that resides on the GPU device
 *
 * This is a fixed-size array that is designed to interact with the Denovo
 * field classes for host<->device transfers.
 *
 * The \c is_initialized() method is useful for error checking. Whenever an
 * instance is constructed by copying another field, or if the \c assign method
 * is called, or the mutable \c data() accessor is called, the instance
 * considers itself initialized.  (The mutable accessor might be used by
 * intializing the data from a kernel call.) If the const \c data() accessor is
 * called before the vector is initialized, a DBC error will result. This is
 * provides some error checking on the host side before passing data to
 * the kernel.
 *
 * \tparam Arch_Switch Whether to use actual GPU or pretend GPU
 * \tparam T           Storage type
 */
//===========================================================================//
template <typename Arch_Switch, typename T>
class Device_Vector
{
  private:
    // Only use specializations!
    Device_Vector();
};

//===========================================================================//
// DEVICE SPECIALIZATION
//===========================================================================//
//! Specialization on Device (GPU storage)
template <typename T>
class Device_Vector<arch::Device, T>
{
    typedef Device_Vector<arch::Device, T> This;
  public:
    typedef arch::Device Arch_t;

    //@{
    //! Container typedefs
    typedef T           value_type;
    typedef T*          pointer;
    typedef const T*    const_pointer;
    typedef std::size_t size_type;
    //@}

    //@{
    //! Field/vector typedefs
    typedef profugus::const_View_Field<T> const_View_Field_t;
    typedef profugus::View_Field<T>       View_Field_t;
    typedef cuda::Host_Vector<T>         Host_Vector_t;
    typedef cuda::Stream<Arch_t>         Stream_t;
    //@}

  public:
    // Construct with number of elements (no initialization!)
    explicit Device_Vector(size_t count);

    // Construct with size and data from a host vector
    explicit Device_Vector(const_View_Field_t hostvec);

    // Construct with size and data from a host vector
    explicit Device_Vector(const Host_Vector_t& hostvec);

    // Copy constructor from another GPU vector (explicit only!)
    explicit Device_Vector(const This& rhs);

    // Free memory on destruct
    ~Device_Vector();

    // Assign same-sized host memory
    void assign(const_View_Field_t hostvec);

    // Assign from host vector (very fast if hostvec is mapped)
    void assign(const Host_Vector_t& hostvec);

    // Assign from another device vector (very fast; on-device)
    void assign(const This& rhs);

    // Asynchronous assignment from host vector
    void assign_async(const Host_Vector_t& hostvec, Stream_t& stream);

    // Swap data with another same-sized GPU vector
    void swap(This& rhs);

    //! Number of stored elements
    size_type size() const { return d_size; }

    //! Whether our data is probably initialized
    bool is_initialized() const { return d_is_initialized; }

    // >>> TRANSFER BACK TO CPU

    // Copy to a host vector
    void to_host(profugus::View_Field<T> out) const;

    // Copy to a host vector
    void to_host(Host_Vector<T>& out) const;

    // Copy to a host vector asynchronously
    void to_host_async(Host_Vector<T>& out, Stream_t& stream) const;

  public:
    /*!
     * \brief Access the underlying data for read/write on the GPU
     *
     * This is needed by kernel call wrappers. It sets the initialized flag.
     *
     * \warning This will always return a pointer to data on the device. Never
     * try to access the underlying data in host (CPU) code. To assign data
     * from the host, use the assign() function.
     */
    pointer data()
    {
        d_is_initialized = true;
        return d_data;
    }

    /*!
     * \brief Access the underlying data for read-only on the GPU
     *
     * \warning This will always return a pointer to data on the device. Never
     * try to access the underlying data in host (CPU) code. To copy the data
     * to the host, use the \c device_to_host function.
     */
    const_pointer cdata() const
    {
        Require(d_is_initialized);
        return d_data;
    }
    const_pointer data() const { return cdata(); }

  private:
    //! Allocation implementation
    void allocate();

  private:
    //! Prevent assignment
    This& operator=(const This& rhs);

  private:
    //! Number of elements
    size_type d_size;

    //! CUDA-allocated data
    pointer d_data;

    //! Whether we've been initialized; for error checking only
    bool d_is_initialized;
};

#ifndef __NVCC__
//===========================================================================//
//! Specialization on Host (CPU storage; GPU emulation)
template <typename T>
class Device_Vector<arch::Host, T>
{
    typedef Device_Vector<arch::Host, T> This;
  public:
    typedef arch::Host Arch_t;

    //@{
    //! Container typedefs
    typedef T           value_type;
    typedef T*          pointer;
    typedef const T*    const_pointer;
    typedef std::size_t size_type;
    //@}

    //@{
    //! Field/vector typedefs
    typedef profugus::const_View_Field<T> const_View_Field_t;
    typedef profugus::View_Field<T>       View_Field_t;
    typedef cuda::Host_Vector<T>         Host_Vector_t;
    typedef cuda::Stream<Arch_t>         Stream_t;
    //@}

  public:
    // Construct with number of elements
    explicit Device_Vector(size_t count);

    // Construct with size and data from a host vector
    explicit Device_Vector(const_View_Field_t hostvec);

    // Construct with size and data from a host vector
    explicit Device_Vector(const Host_Vector_t& hostvec);

    // Copy constructor from another GPU vector (explicit only!)
    explicit Device_Vector(const This& rhs);

    // Assign same-sized host memory
    void assign(const_View_Field_t hostvec);

    // Assign from host vector (very fast if hostvec is mapped)
    void assign(const Host_Vector_t& hostvec);

    // Copy from "device"
    void assign(const This& rhs);

    // "Asynchronous" assignment from host vector
    void assign_async(const Host_Vector_t& hostvec, Stream_t&)
    {
        this->assign(hostvec);
    }

    // Swap data with another same-sized GPU vector
    void swap(This& rhs);

    //! Number of stored elements
    size_type size() const { return d_storage.size(); }

    //! Whether our data is probably initialized
    bool is_initialized() const { return d_is_initialized; }

    // >>> TRANSFER BACK TO CPU

    // Copy to a host vector
    void to_host(profugus::View_Field<T> out) const;

    // Copy to a host vector
    void to_host(Host_Vector<T>& out) const;

    // "Asynchronous" copy to host vector
    void to_host_async(Host_Vector_t& out, Stream_t&) const
    {
        this->to_host(out);
    }

  public:
    /*!
     * \brief Access the underlying data for read/write in fake GPU kernels
     *
     * This is needed by kernel call wrappers. It sets the initialized flag.
     */
    pointer data()
    {
        d_is_initialized = true;
        return &d_storage[0];
    }

    /*!
     * \brief Access the underlying data for read-only in fake GPU kernels
     */
    const_pointer cdata() const
    {
        Require(d_is_initialized);
        return &d_storage[0];
    }
    const_pointer data() const { return cdata(); }

  private:
    //! Prevent assignment
    This& operator=(const This& rhs);

  private:
    //! Underlying storage container
    std::vector<T> d_storage;

    //! Whether we've been initialized; for error checking only
    bool d_is_initialized;
};
#endif

//===========================================================================//
template <typename Arch_Switch, typename T>
inline void device_to_host(
        const Device_Vector<Arch_Switch, T>& in,
        profugus::View_Field<T> out)
{
    in.to_host(out);
}

//---------------------------------------------------------------------------//
template <typename Arch_Switch, typename T>
inline void device_to_host(
        const Device_Vector<Arch_Switch, T>& in,
        Host_Vector<T>& out)
{
    in.to_host(out);
}

//---------------------------------------------------------------------------//
template <typename Arch_Switch, typename T>
inline void device_to_host_async(
        const Device_Vector<Arch_Switch,T>& in,
        Host_Vector<T>&                     out,
        Stream<Arch_Switch>&                stream)
{
    in.to_host_async(out, stream);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Swap two same-sized vectors on the device by switching pointers.
 *
 * This is fast as it performs no copy operations.
 */
template <typename Arch_Switch, typename T>
inline void swap(
        Device_Vector<Arch_Switch, T>& left,
        Device_Vector<Arch_Switch, T>& right)
{
    left.swap(right);
}

//===========================================================================//

} // end namespace cuda

#endif // cuda_utils_Device_Vector_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Device_Vector.hh
//---------------------------------------------------------------------------//
