// vim: set ft=cuda: ---------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Texture_Vector_host.t.hh
 * \author Seth R Johnson
 * \date   Fri Sep 20 10:08:43 2013
 * \brief  Texture_Vector template member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Texture_Vector_host_t_hh
#define cuda_utils_Texture_Vector_host_t_hh

#include "Texture_Vector.hh"

#include "harness/DBC.hh"
#include "utils/View_Field.hh"
#include "Host_Vector.hh"
#include "Device_Vector.hh"

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
Texture_Vector<arch::Host,T>::Texture_Vector(size_t count)
  : d_data(count)
  , d_is_initialized(false)
{
    REQUIRE(count > 0);

    ENSURE(size() == count);
    ENSURE(!is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct by copying data from the host in pageable memory
 */
template<typename T>
Texture_Vector<arch::Host,T>::Texture_Vector(const_View_Field_t hostvec)
  : d_data(hostvec.begin(), hostvec.end())
  , d_is_initialized(true)
{
    ENSURE(size() == hostvec.size());
    ENSURE(is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct by copying data from the host in pageable memory
 */
template<typename T>
Texture_Vector<arch::Host,T>::Texture_Vector(const Host_Vector_t& hostvec)
  : d_data(hostvec.begin(), hostvec.end())
  , d_is_initialized(true)
{
    ENSURE(size() == hostvec.size());
    ENSURE(is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct with data stored in a GPU vector.
 */
template<typename T>
Texture_Vector<arch::Host,T>::Texture_Vector(const Device_Vector_t& devicevecec)
  : d_data(devicevecec.data(), devicevecec.data() + devicevecec.size())
  , d_is_initialized(true)
{
    ENSURE(size() == devicevecec.size());
    ENSURE(is_initialized());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Clean up texture on destruct
 */
template<typename T>
Texture_Vector<arch::Host,T>::~Texture_Vector()
{
    /* * */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign CUDA texture object
 */
template<typename T>
void Texture_Vector<arch::Host,T>::assign(const_View_Field_t hostvec)
{
    REQUIRE(hostvec.size() == size());
    d_data.assign(hostvec.begin(), hostvec.end());
    d_is_initialized = true;
    ENSURE(is_initialized());
}

//---------------------------------------------------------------------------//
template<typename T>
void Texture_Vector<arch::Host,T>::assign(const Host_Vector_t& hostvec)
{
    REQUIRE(hostvec.size() == size());
    d_data.assign(hostvec.begin(), hostvec.end());
    d_is_initialized = true;
    ENSURE(is_initialized());
}

//---------------------------------------------------------------------------//
template<typename T>
void Texture_Vector<arch::Host,T>::assign(const Device_Vector_t& devicevec)
{
    REQUIRE(devicevec.size() == size());
    d_data.assign(devicevec.data(), devicevec.data() + devicevec.size());
    d_is_initialized = true;
    ENSURE(is_initialized());
}
//---------------------------------------------------------------------------//
} // end namespace cuda

#endif // cuda_utils_Texture_Vector_host_t_hh

//---------------------------------------------------------------------------//
//                 end of Texture_Vector_host.t.hh
//---------------------------------------------------------------------------//
