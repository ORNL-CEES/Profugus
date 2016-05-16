//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Multi_Vector.t.hh
 * \author Seth R Johnson
 * \date   Fri Aug 02 10:24:45 2013
 * \brief  Multi_Vector inline member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Multi_Vector_t_hh
#define CudaUtils_cuda_utils_Multi_Vector_t_hh

#include "Multi_Vector.hh"
#include "utils/View_Field.hh"

namespace cuda
{
//---------------------------------------------------------------------------//
template <typename Arch_Switch, typename T>
Multi_Vector<Arch_Switch,T>::Multi_Vector(size_type count)
  : d_device_vectors(count, NULL)
{
    ENSURE(size() == count);
}
//---------------------------------------------------------------------------//
/*!
 * \brief Clean up on destruction
 */
template <typename Arch_Switch, typename T>
Multi_Vector<Arch_Switch,T>::~Multi_Vector()
{
    for (typename Vec_P_DV::iterator it = d_device_vectors.begin(),
            end_it = d_device_vectors.end();
            it != end_it;
            ++it)
    {
        Device_Vector_t* dv = *it;
        if (dv != NULL)
            delete dv;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assign data to an index for the first time
 */
template <typename Arch_Switch, typename T>
void Multi_Vector<Arch_Switch,T>::initialize(
        size_type index,
        const_View_Field_t data)
{
    REQUIRE(index < size());
    REQUIRE(!is_initialized(index));
    REQUIRE(data.size() > 0);

    d_device_vectors[index] = new Device_Vector_t(data);
    CHECK(d_device_vectors[index]->is_initialized());

    ENSURE(is_initialized(index));
}

//---------------------------------------------------------------------------//
} // end namespace cuda

#endif // CudaUtils_cuda_utils_Multi_Vector_t_hh

//---------------------------------------------------------------------------//
//                   end of cuda_utils/Multi_Vector.t.hh
//---------------------------------------------------------------------------//
