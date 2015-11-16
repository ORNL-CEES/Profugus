//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_utils/Texture_Vector_Kernel.cuh
 * \author Seth R Johnson
 * \date   Fri Sep 20 18:06:54 2013
 * \brief  Texture_Vector_Kernel kernel declarations.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Texture_Vector_Kernel_cuh
#define cuda_utils_Texture_Vector_Kernel_cuh

#include "Utils/harness/DBC.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
// HOST INTERFACES
//---------------------------------------------------------------------------//

// Declare class so it can be our friend
template<typename Arch_Switch, typename T> class Texture_Vector;

//===========================================================================//
// CLASS DEFINITIONS
//===========================================================================//
/*!
 * \class Texture_Vector_Kernel
 * \brief Emulate vector accesses using texture memory
 *
 * The Texture_Vector object must be created on the host, and then passed
 * to the kernel.
 *
 * This does not do any size checking on accesses.
 */
//===========================================================================//
template <typename Arch_Switch, typename T>
class Texture_Vector_Kernel
{
  private:
    // Only use specializations!
    Texture_Vector_Kernel();
};

#if defined(USE_CUDA) && !defined(PSEUDO_CUDA)
//===========================================================================//
//! Specialization on Device (GPU implementation)
template <typename T>
class Texture_Vector_Kernel<arch::Device, T>
{
    typedef Texture_Vector_Kernel<arch::Device, T> This;
  public:
    typedef arch::Device Arch_t;

    typedef T      value_type;
    typedef size_t size_type;
  private:
    // >>> DATA

    //! Texture object
    cudaTextureObject_t d_texture;

  private:
    // >>> CONSTRUCTION

    friend class Texture_Vector<Arch_t, value_type>;

    // Initialize on the host (only our friend can create)
    Texture_Vector_Kernel(cudaTextureObject_t tex) : d_texture(tex)
    {
        /* * */
    }

  public:
    //! Implicit copy constructor during kernel call
    Texture_Vector_Kernel(const This& rhs) : d_texture(rhs.d_texture)
    {
        /* * */
    }

#ifdef __NVCC__
    //! Access texture memory like a vector
    __device__ __inline__ value_type operator[] (int i) const
    {
        return tex1Dfetch<value_type>(d_texture, i);
    }
#endif
};

#ifdef __NVCC__
//! Specialization on double
template <>
__device__ __inline__ double
Texture_Vector_Kernel<arch::Device,double>::operator[](int i) const
{
    int2 v = tex1Dfetch<int2>(d_texture, i);
    return __hiloint2double(v.y, v.x);
}
#endif

#endif // USE_CUDA && !PSEUDO_CUDA

//===========================================================================//
//! Specialization for CPU emulation code
template <typename T>
class Texture_Vector_Kernel<arch::Host, T>
{
    typedef Texture_Vector_Kernel<arch::Host, T> This;
  public:
    typedef arch::Host Arch_t;

    typedef T      value_type;
    typedef size_t size_type;

  private:
    // >>> DATA

    //! Pointer to underlying data
    const value_type* d_data;

    //! Size for error checking
    size_type d_size;

  private:
    // >>> CONSTRUCTION

    friend class Texture_Vector<Arch_t, value_type>;

    // Initialize on the host (only our friend can create)
    Texture_Vector_Kernel(
            const value_type* data,
            size_type   size)
      : d_data(data)
      , d_size(size)
    {
        Require(d_data);
    }

  public:
    //! Implicit copy constructor during kernel call
    Texture_Vector_Kernel(const This& rhs)
      : d_data(rhs.d_data)
      , d_size(rhs.d_size)
    {
        /* * */
    }

    //! Access texture memory like a vector
    value_type operator[] (int i) const
    {
        Require(i >= 0 && i < d_size);
        return d_data[i];
    }
};
//---------------------------------------------------------------------------//
} // end namespace cuda
#endif // cuda_utils_Texture_Vector_Kernel_cuh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/Texture_Vector_Kernel.cuh
//---------------------------------------------------------------------------//
