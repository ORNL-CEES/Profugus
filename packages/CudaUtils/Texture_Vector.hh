//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Texture_Vector.hh
 * \author Seth R Johnson
 * \date   Fri Sep 20 10:08:43 2013
 * \brief  Texture_Vector class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Texture_Vector_hh
#define cuda_utils_Texture_Vector_hh

#include <config.h>
#include <vector>
#if defined(USE_CUDA) && !defined(PSEUDO_CUDA)
#include <cuda_runtime.h>
#endif // USE_CUDA && !PSEUDO_CUDA

#include "Definitions.hh"
#include "Texture_Vector_Kernel.cuh"

// Declare fields to avoid header propagation into CUDA kernel files
namespace profugus
{
template<typename T> class View_Field;
template<typename T> class const_View_Field;
}

namespace cuda
{
// Declare vector classes
template <typename Arch_Switch, typename T> class Device_Vector;
template <typename T> class Host_Vector;

//===========================================================================//
/*!
 * \class Texture_Vector
 * \brief A vector-like interface for accessing read-only data through textures
 *
 * Textures are accessed and cached using special hardware on the GPU. They can
 * only be assigned once, but they can be swapped (fast, constant-time
 * operation) with other texture vector objects.
 */
//===========================================================================//
template <typename Arch_Switch, typename T>
class Texture_Vector
{
  private:
    // Only use specializations!
    Texture_Vector();
};

#if defined(USE_CUDA) && !defined(PSEUDO_CUDA)
//===========================================================================//
// DEVICE SPECIALIZATION
//===========================================================================//
//! Specialization on Device (GPU storage)
template <typename T>
class Texture_Vector<arch::Device, T>
{
    typedef Texture_Vector<arch::Device, T> This;

  public:
    typedef arch::Device Arch_t;

    //@{
    //! Container typedefs
    typedef T        value_type;
    typedef T*       pointer;
    typedef const T* const_pointer;
    typedef size_t   size_type;
    //@}

    //@{
    //! Field/vector typedefs
    typedef profugus::const_View_Field<T>  const_View_Field_t;
    typedef profugus::View_Field<T>        View_Field_t;
    typedef cuda::Host_Vector<T>          Host_Vector_t;
    typedef cuda::Device_Vector<Arch_t,T> Device_Vector_t;
    //@}

    typedef Texture_Vector_Kernel<Arch_t, value_type> TVK_t;

  public:
    // >>> CONSTRUCTION

    // Construct with number of elements (no initialization!)
    explicit Texture_Vector(size_t count);

    // Construct with size and data from a host vector
    explicit Texture_Vector(const_View_Field_t hostvec);

    // Construct with size and data from a host vector
    explicit Texture_Vector(const Host_Vector_t& hostvec);

    // Construct with data already on the device
    explicit Texture_Vector(const Device_Vector_t& devicevec);

    // Clean up texture
    ~Texture_Vector();

    // >>> ACCESSING ON HOST

    // Assign same-sized host memory
    void assign(const_View_Field_t hostvec);

    // Assign from host vector
    void assign(const Host_Vector_t& hostvec);

    // Assign from device vector
    void assign(const Device_Vector_t& devicevec);

    //! Number of stored elements
    size_type size() const { return d_data.size(); }

    //! Whether our data is probably initialized
    bool is_initialized() const { return d_is_initialized; }

    // >>> ACCESSING ON DEVICE

    // Return the kernel accessor by value
    TVK_t data() const
    {
        return TVK_t(d_texture);
    }

  private:
    // Prevent copying!
    Texture_Vector(const This& rhs);
    This& operator=(const This& rhs);

  private:
    // >>> IMPLEMENTATION

    //! Object-agnostic assigment
    template <class U> void assign_any(const U& vec);

    //! Create the texture object
    void create_texture();

    //! Destroy the texture object
    void destroy_texture();

  private:
    // >>> IMPLEMENTATION

    //! Storage for data associated with texture
    Device_Vector_t d_data;

    //! Underlying texture object
    cudaTextureObject_t d_texture;

    //! Whether we've been initialized (for error checking)
    bool d_is_initialized;
};
#endif // USE_CUDA && !PSEUDO_CUDA

//===========================================================================//
//! Specialization on Device (CPU storage)
template <typename T>
class Texture_Vector<arch::Host, T>
{
    typedef Texture_Vector<arch::Host, T> This;

  public:
    typedef arch::Host Arch_t;

    //@{
    //! Container typedefs
    typedef T        value_type;
    typedef T*       pointer;
    typedef const T* const_pointer;
    typedef size_t   size_type;
    //@}

    //@{
    //! Field/vector typedefs
    typedef profugus::const_View_Field<T>  const_View_Field_t;
    typedef profugus::View_Field<T>        View_Field_t;
    typedef cuda::Host_Vector<T>          Host_Vector_t;
    typedef cuda::Device_Vector<Arch_t,T> Device_Vector_t;
    //@}

    typedef Texture_Vector_Kernel<Arch_t, value_type> TVK_t;

  public:
    // >>> CONSTRUCTION

    // Construct with number of elements (no initialization!)
    explicit Texture_Vector(size_t count);

    // Construct with size and data from a host vector
    explicit Texture_Vector(const_View_Field_t hostvec);

    // Construct with size and data from a host vector
    explicit Texture_Vector(const Host_Vector_t& hostvec);

    // Construct with data already on the device
    explicit Texture_Vector(const Device_Vector_t& devicevec);

    // Clean up texture
    ~Texture_Vector();

    // >>> ACCESSING ON HOST

    // Assign same-sized host memory
    void assign(const_View_Field_t hostvec);

    // Assign from host vector
    void assign(const Host_Vector_t& hostvec);

    // Assign from device vector
    void assign(const Device_Vector_t& devicevec);

    //! Number of stored elements
    size_type size() const { return d_data.size(); }

    //! Whether our data is probably initialized
    bool is_initialized() const { return d_is_initialized; }

    // >>> ACCESSING ON DEVICE

    // Return the kernel accessor by value
    TVK_t data() const
    {
        Require(is_initialized());

        return TVK_t(&d_data.front(), d_data.size());
    }

  private:
    // Prevent copying
    Texture_Vector(const This& rhs);
    This& operator=(const This& rhs);

  private:
    // >>> IMPLEMENTATION

    //! Storage for data
    std::vector<value_type> d_data;

    //! Whether we've been initialized
    bool d_is_initialized;
};

//===========================================================================//
/*!
 * \brief Swap two same-sized textures on the device by switching pointers.
 *
 * This is fast as it performs no copy operations.
 */
template <typename Arch_Switch, typename T>
inline void swap(
        Texture_Vector<Arch_Switch, T>& left,
        Texture_Vector<Arch_Switch, T>& right)
{
    left.swap(right);
}

//===========================================================================//
} // end namespace cuda

#endif // cuda_utils_Texture_Vector_hh

//---------------------------------------------------------------------------//
//                 end of Texture_Vector.hh
//---------------------------------------------------------------------------//
