//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Device_Vector_host.t.hh
 * \author Seth R Johnson
 * \date   Thu Aug 01 11:33:12 2013
 * \brief  Device_Vector template member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Device_Vector_host_t_hh
#define cuda_utils_Device_Vector_host_t_hh

#include "Device_Vector.hh"

#include <cstdlib>
#include <utility>

#include "harness/DBC.hh"
#include "utils/View_Field.hh"
#include "Host_Vector.hh"

namespace cuda_utils
{
//---------------------------------------------------------------------------//
/*!
 * \brief Construct with number of elements
 *
 * This does *not* perform any initialization, which might require a kernel
 * call.
 */
template<typename T>
Device_Vector<arch::Host,T>::Device_Vector(size_t count)
  : d_storage(count)
  , d_is_initialized(false)
{
    REQUIRE(count > 0);
    ENSURE(!is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct by copying data from the host in pageable memory
 */
template<typename T>
Device_Vector<arch::Host,T>::Device_Vector(const_View_Field_t hostvec)
  : d_storage(hostvec.begin(), hostvec.end())
  , d_is_initialized(true)
{
    REQUIRE(hostvec.size() > 0);
    ENSURE(is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct by copying data from the host in pageable memory
 */
template<typename T>
Device_Vector<arch::Host,T>::Device_Vector(const Host_Vector_t& hostvec)
  : d_storage(hostvec.begin(), hostvec.end())
  , d_is_initialized(true)
{
    REQUIRE(hostvec.size() > 0);
    ENSURE(is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct with a copy of data in another GPU vector.
 *
 * This is an \c explicit constructor because we don't want to be unwittingly
 * doing memory copy operations on the GPU.
 *
 * \note If DBC is on, an assertion will be thrown if the GPU vector being
 * copied is uninitialized.
 */
template<typename T>
Device_Vector<arch::Host,T>::Device_Vector(const This& rhs)
  : d_storage(rhs.d_storage)
  , d_is_initialized(rhs.d_is_initialized)
{
    REQUIRE(rhs.is_initialized());
    ENSURE(is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign same-sized host memory
 */
template<typename T>
void Device_Vector<arch::Host,T>::assign(const_View_Field_t hostvec)
{
    REQUIRE(hostvec.size() == size());
    d_storage.assign(hostvec.begin(), hostvec.end());
    d_is_initialized = true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign same-sized host memory from a Host_Vector
 */
template<typename T>
void Device_Vector<arch::Host,T>::assign(const Host_Vector_t& hostvec)
{
    REQUIRE(hostvec.size() == size());
    d_storage.assign(hostvec.begin(), hostvec.end());
    d_is_initialized = true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign same-sized host memory from a Host_Vector
 */
template<typename T>
void Device_Vector<arch::Host,T>::assign(const This& rhs)
{
    REQUIRE(rhs.size() == size());
    d_storage.assign(rhs.d_storage.begin(), rhs.d_storage.end());
    d_is_initialized = true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Swap two same-sized device vectors
 */
template<typename T>
void Device_Vector<arch::Host,T>::swap(This& rhs)
{
    std::swap(d_storage, rhs.d_storage);
    std::swap(d_is_initialized, rhs.d_is_initialized);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do a "device"-to-host transfer for standard memory
 */
template<typename T>
void Device_Vector<arch::Host,T>::to_host(profugus::View_Field<T> out) const
{
    REQUIRE(size() == out.size());
    REQUIRE(is_initialized());

    std::memcpy(&out[0], data(), size() * sizeof(T));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do a "device"-to-host transfer for pinned memory
 */
template<typename T>
void Device_Vector<arch::Host,T>::to_host(Host_Vector<T>& out) const
{
    REQUIRE(size() == out.size());
    REQUIRE(is_initialized());

    std::memcpy(&out[0], data(), size() * sizeof(T));
}

//===========================================================================//

} // end namespace cuda_utils

#endif // cuda_utils_Device_Vector_host_t_hh

//---------------------------------------------------------------------------//
//                   end of cuda_utils/Device_Vector_host.t.hh
//---------------------------------------------------------------------------//
