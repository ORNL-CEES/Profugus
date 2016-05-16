//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/Multi_Vector.hh
 * \author Seth R Johnson
 * \date   Fri Aug 02 10:24:45 2013
 * \brief  Multi_Vector class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_cuda_utils_Multi_Vector_hh
#define CudaUtils_cuda_utils_Multi_Vector_hh

#include <vector>
#include "harness/DBC.hh"
#include "Device_Vector.hh"

namespace cuda
{

//===========================================================================//
/*!
 * \class Multi_Vector
 * \brief Independent storage of multiple associated device vectors
 *
 * This is used, for example, when caching group-dependent ray segment optical
 * thicknesses in MOC. Not all entries have to be assigned, and the vectors are
 * almost certainly not going to be contiguous on the GPU.
 *
 * The class's utility is primarily the automated memory management: there's no
 * need to manually free a bunch of device vectors.
 *
 * \tparam Vector_T A Device_Vector-like object.
 */
/*!
 * \example cuda_utils/test/tstMulti_Vector.cc
 *
 * Test of Multi_Vector.
 */
//===========================================================================//

template <typename Arch_Switch, typename T>
class Multi_Vector
{
    typedef Multi_Vector<Arch_Switch,T> This;
  public:
    typedef Arch_Switch Arch_t;

    //@{
    //! Container typedefs
    typedef T           value_type;
    typedef std::size_t size_type;
    //@}

    typedef Device_Vector<Arch_t, value_type>     Device_Vector_t;
    typedef profugus::const_View_Field<value_type> const_View_Field_t;
    typedef std::vector<Device_Vector_t*>         Vec_P_DV;

  private:
    //! Underlying storage
    Vec_P_DV d_device_vectors;

  public:
    /*!
     * \brief Construct an MDV with the number of vectors it can hold
     *
     * This is a fairly lightweight operation: it allocates sizeof(P) * count
     * memory.
     */
    explicit Multi_Vector(size_type count);

    // Clean up on destruction
    ~Multi_Vector();

    //! Length of multivector (highest index + 1)
    size_type size() const { return d_device_vectors.size(); }

    //! Whether an index has been assigned
    bool is_initialized(size_type index) const
    {
        REQUIRE(index < size());
        return d_device_vectors[index] != NULL;
    }

    // Assign data to an index for the first time
    void initialize(size_type index, const_View_Field_t data);

    // Return device vector for an index; raise DBC error if not yet assigned
    const Device_Vector_t& operator[](size_type index) const
    {
        REQUIRE(index < size());
        REQUIRE(is_initialized(index));
        return *d_device_vectors[index];
    }

    Device_Vector_t& operator[](size_type index)
    {
        REQUIRE(index < size());
        REQUIRE(is_initialized(index));
        return *d_device_vectors[index];
    }

  private:
    // Disallow copying for now
    Multi_Vector(const Multi_Vector&);
};

} // end namespace cuda

//---------------------------------------------------------------------------//
#endif // CudaUtils_cuda_utils_Multi_Vector_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Multi_Vector.hh
//---------------------------------------------------------------------------//
