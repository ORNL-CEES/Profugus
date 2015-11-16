// vim: set ft=cuda: ---------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Texture_Vector_device.t.hh
 * \author Seth R Johnson
 * \date   Fri Sep 20 10:08:43 2013
 * \brief  Texture_Vector template member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Texture_Vector_device_t_hh
#define cuda_utils_Texture_Vector_device_t_hh

#include "Texture_Vector.hh"

#include <cuda_runtime.h>
#include <cstdlib>

#include "Utils/harness/DBC.hh"
#include "Utils/comm/Logger.hh"
#include "Utils/utils/View_Field.hh"
#include "Host_Vector.hh"
#include "Device_Vector.hh"

#include "CudaDBC.hh"

//---------------------------------------------------------------------------//
// ANONYMOUS HELPER FUNCTIONS
//---------------------------------------------------------------------------//
namespace
{
template<typename T>
struct Resource_Desc_Data;

template<>
struct Resource_Desc_Data<float>
{
    static cudaChannelFormatKind kind() { return cudaChannelFormatKindFloat; }
    static int x_bytes() { return 4; }
    static int y_bytes() { return 0; }
};

template<>
struct Resource_Desc_Data<int>
{
    static cudaChannelFormatKind kind() { return cudaChannelFormatKindSigned; }
    static int x_bytes() { return sizeof(int); }
    static int y_bytes() { return 0; }
};

template<>
struct Resource_Desc_Data<unsigned int>
{
    static cudaChannelFormatKind kind() { return cudaChannelFormatKindUnsigned; }
    static int x_bytes() { return sizeof(unsigned int); }
    static int y_bytes() { return 0; }
};

template<>
struct Resource_Desc_Data<double>
{
    // Double gets stored as two integer elements
    static cudaChannelFormatKind kind() { return cudaChannelFormatKindSigned; }
    static int x_bytes() { return 4; }
    static int y_bytes() { return 4; }
};

}
//---------------------------------------------------------------------------//

namespace cuda
{
//---------------------------------------------------------------------------//
// PUBLIC METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Construct with number of elements
 *
 * This does *not* perform any initialization.
 */
template<typename T>
Texture_Vector<arch::Device,T>::Texture_Vector(size_t count)
  : d_data(count)
  , d_is_initialized(false)
{
    Require(count > 0);

    Ensure(size() == count);
    Ensure(!is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct by copying data from the host in pageable memory
 */
template<typename T>
Texture_Vector<arch::Device,T>::Texture_Vector(const_View_Field_t hostvec)
  : d_data(hostvec)
  , d_is_initialized(false)
{
    this->create_texture();

    Ensure(size() == hostvec.size());
    Ensure(is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct by copying data from the host in pageable memory
 */
template<typename T>
Texture_Vector<arch::Device,T>::Texture_Vector(const Host_Vector_t& hostvec)
  : d_data(hostvec)
  , d_is_initialized(false)
{
    this->create_texture();

    Ensure(size() == hostvec.size());
    Ensure(is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct with data stored in a GPU vector.
 */
template<typename T>
Texture_Vector<arch::Device,T>::Texture_Vector(const Device_Vector_t& devvec)
  : d_data(devvec)
  , d_is_initialized(false)
{
    this->create_texture();

    Ensure(size() == devvec.size());
    Ensure(is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Clean up texture on destruct
 *
 * The GPU memory is automatically freed by the Device_Vector destructor.
 */
template<typename T>
Texture_Vector<arch::Device,T>::~Texture_Vector()
{
    if (is_initialized())
    {
        try
        {
            this->destroy_texture();
        }
        catch (const profugus::assertion& e)
        {
            log(profugus::WARNING)
                << "Error: failed to destroy texture: " << e.what();
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign CUDA texture object
 *
 * This is more expensive than copying memory because it requires
 *
 * \warning Don't modify a texture while the kernel is in use.
 */
template<typename T>
void Texture_Vector<arch::Device,T>::assign(const_View_Field_t hostvec)
{
    Require(hostvec.size() == size());
    this->assign_any(hostvec);
    Ensure(is_initialized());
}

//---------------------------------------------------------------------------//
template<typename T>
void Texture_Vector<arch::Device,T>::assign(const Host_Vector_t& hostvec)
{
    Require(hostvec.size() == size());
    this->assign_any(hostvec);
    Ensure(is_initialized());
}

//---------------------------------------------------------------------------//
template<typename T>
void Texture_Vector<arch::Device,T>::assign(const Device_Vector_t& devicevec)
{
    Require(devicevec.size() == size());
    this->assign_any(devicevec);
    Ensure(is_initialized());
}

//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
template<typename T>
template<typename U>
void Texture_Vector<arch::Device,T>::assign_any(const U& vec)
{
    Require(vec.size() == size());

    // If we're already bound to a texture, delete the texture object
    if (is_initialized())
        destroy_texture();

    // Copy data to device (or in device)
    d_data.assign(vec);

    // Initialize texture from device memory
    create_texture();

    Ensure(is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize the CUDA texture object
 */
template<typename T>
void Texture_Vector<arch::Device,T>::create_texture()
{
    Require(d_data.is_initialized());
    Require(!is_initialized());

    // Resource descriptor
    cudaResourceDesc res_desc;
    std::memset(&res_desc, 0, sizeof(res_desc));
    res_desc.resType = cudaResourceTypeLinear;
    res_desc.res.linear.devPtr = d_data.data();

    // Type and number of bits belonging to the elements in the vector
    // (i.e. int is signed, has one element (8 * 4 bits) ;
    //  uint2 is unsigned, has two elements (each 8 * 4 bits)
    typedef Resource_Desc_Data<T> RDD_t;
    res_desc.res.linear.desc.f = RDD_t::kind();
    res_desc.res.linear.desc.x = 8 * RDD_t::x_bytes();
    res_desc.res.linear.desc.y = 8 * RDD_t::y_bytes();
    res_desc.res.linear.sizeInBytes = d_data.size() * sizeof(value_type);

    // Texture descriptor
    cudaTextureDesc tex_desc;
    std::memset(&tex_desc, 0, sizeof(tex_desc));
    tex_desc.readMode = cudaReadModeElementType;

    CudaCall(cudaCreateTextureObject(&d_texture, &res_desc, &tex_desc, NULL));

    d_is_initialized = true;

    Ensure(is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destroy the CUDA texture object
 *
 * This unbinds the texture from the GPU but does not release the device
 * memory. That only happens on destruction.
 *
 * \warning Don't do this while the kernel using this texture is being executed!
 */
template<typename T>
void Texture_Vector<arch::Device,T>::destroy_texture()
{
    Require(is_initialized());

    CudaCall(cudaDestroyTextureObject(d_texture));

    d_is_initialized = false;

    Ensure(!is_initialized());
}

//---------------------------------------------------------------------------//
} // end namespace cuda

#endif // cuda_utils_Texture_Vector_device_t_hh

//---------------------------------------------------------------------------//
//                 end of Texture_Vector_device.t.hh
//---------------------------------------------------------------------------//
