//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Metaclass.t.hh
 * \author Seth R Johnson
 * \date   Fri Aug 30 20:59:34 2013
 * \brief  Metaclass template member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Metaclass_t_hh
#define utils_Metaclass_t_hh

#include "Metaclass.hh"

#include "harness/Warnings.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// STATIC DATA
//---------------------------------------------------------------------------//

//! Metadata about extra entries
template<class I>
typename Metaclass<I>::Vec_Member_Md Metaclass<I>::d_member_md;

//! Accumulated data size (bytes)
template<class I>
unsigned int Metaclass<I>::d_storage_size = 0;

//! Total number of extant instances (for error checking)
template<class I>
unsigned int Metaclass<I>::d_num_instances = 0;

//---------------------------------------------------------------------------//
/*!
 * \brief Access name of a member given a member ID
 */
template<typename I>
std::string Metaclass<I>::name(unsigned int member_id)
{
    Require(member_id < size());

    return d_member_md[member_id].name;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reset the metadata about extra entries
 *
 * \warning This can only be called when no instances exist: otherwise,
 * destructors won't be called and things will be sad.
 */
template<class I>
void Metaclass<I>::reset()
{
#if UTILS_DBC > 0
    Insist(Metaclass::d_num_instances == 0,
            "Reset can only be called when no instances exist.");
#endif

    // Delete all metadata
    Metaclass::d_member_md.clear();
    Metaclass::d_storage_size = 0;
}

//---------------------------------------------------------------------------//
// CLASS METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize, allocating storage
 */
template<class I>
Metaclass<I>::Metaclass()
  : d_data(NULL)
{
#if UTILS_DBC > 0
    // Increment the particle count for error checking
    ++Metaclass::d_num_instances;
#endif

    // Allocate
    this->alloc();

    // Apply constructors
    for (typename Vec_Member_Md::const_iterator md = d_member_md.begin(),
            end_md = d_member_md.end();
            md != end_md;
            ++md)
    {
        Check(md->member_manager);
        char* data = d_data + md->data_offset;
        md->member_manager->construct(reinterpret_cast<void*>(data));
    }

    Ensure(d_data);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor
 *
 * This creates copies of the data as defined by the member managers.
 */
template<class I>
Metaclass<I>::Metaclass(const This& rhs)
  : d_data(NULL)
{
#if UTILS_DBC > 0
    // Increment the particle count for error checking
    ++Metaclass::d_num_instances;
#endif
    // Allocate
    this->alloc();

    // Apply constructors
    for (typename Vec_Member_Md::const_iterator md = d_member_md.begin(),
            end_md = d_member_md.end();
            md != end_md;
            ++md)
    {
        Check(md->member_manager);
        char*           data = d_data + md->data_offset;
        const char* rhs_data = rhs.d_data + md->data_offset;
        md->member_manager->copy_construct(
                reinterpret_cast<void*>(data),
                reinterpret_cast<const void*>(rhs_data));
    }

    Ensure(d_data);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor
 */
template<class I>
Metaclass<I>::Metaclass(const char* const buffer_begin, unsigned int size)
  : d_data(NULL)
{
#if UTILS_DBC > 0
    // Increment the particle count for error checking
    ++Metaclass::d_num_instances;
#endif
    // Allocate
    this->alloc();

    // Beginning location for the buffer to unpack an object
    const char* buffer = buffer_begin;

    // Apply constructors
    for (typename Vec_Member_Md::const_iterator md = d_member_md.begin(),
            end_md = d_member_md.end();
            md != end_md;
            ++md)
    {
        Check(md->member_manager);
        char* data = d_data + md->data_offset;
        buffer = md->member_manager->unpack_construct(
                reinterpret_cast<void*>(data),
                buffer);
        Check(buffer <= buffer_begin + size);
    }
    Check(buffer == buffer_begin + size);
    Ensure(d_data);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment
 *
 * This creates *shallow* copies of the data.
 */
template<class I>
Metaclass<I>& Metaclass<I>::operator=(const This& rhs)
{
    Require(d_data);
    Require(rhs.d_data);

    // Apply assignment
    for (typename Vec_Member_Md::const_iterator md = d_member_md.begin(),
            end_md = d_member_md.end();
            md != end_md;
            ++md)
    {
        Check(md->member_manager);
        char*           data = d_data + md->data_offset;
        const char* rhs_data = rhs.d_data + md->data_offset;
        md->member_manager->assign(
                reinterpret_cast<void*>(data),
                reinterpret_cast<const void*>(rhs_data));
    }

    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Clean up extra data on destruction
 *
 * This calls the "destructor" on each piece of data in there.
 */
template<class I>
Metaclass<I>::~Metaclass()
{
    for (typename Vec_Member_Md::const_iterator md = d_member_md.begin(),
            end_md = d_member_md.end();
            md != end_md;
            ++md)
    {
        try
        {
            Check(md->member_manager);
            char* data = d_data + md->data_offset;
            md->member_manager->destroy(reinterpret_cast<void*>(data));
        }
        catch (const profugus::assertion& err)
        {
            ADD_WARNING("Caught exception while deallocating entry '"
                    << md->name << "'");
        }
    }

    // Free memory
    std::free(d_data);

#if UTILS_DBC > 0
    // Decrement the extant instance count
    --Metaclass::d_num_instances;
#endif
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the packed size of this metaclass.
 *
 * This calls the packed size function for each element of data.
 */
template<class I>
typename Metaclass<I>::size_type Metaclass<I>::packed_size() const
{
    size_type result = 0;

    for (typename Vec_Member_Md::const_iterator md = d_member_md.begin(),
            end_md = d_member_md.end();
            md != end_md;
            ++md)
    {
        Check(md->member_manager);
        char* data = d_data + md->data_offset;
        // Add the packed size of this data object
        result += md->member_manager->packed_size(
                reinterpret_cast<void*>(data));
    }

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack into a properly-sized buffer
 *
 * \return Pointer to the new buffer location
 */
template<class I>
char* Metaclass<I>::pack(char* buffer) const
{
    for (typename Vec_Member_Md::const_iterator md = d_member_md.begin(),
            end_md = d_member_md.end();
            md != end_md;
            ++md)
    {
        Check(md->member_manager);
        char* data = d_data + md->data_offset;
        // Pack into the buffer and set the newly incremented pointer
        buffer = md->member_manager->pack(
                reinterpret_cast<void*>(data),
                buffer);
    }
    return buffer;
}

//---------------------------------------------------------------------------//
// PRIVATE METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Allocate memory
 */
template<typename I>
void Metaclass<I>::alloc()
{
    d_data = reinterpret_cast<char*>(std::malloc(storage_size()));

    Ensure(d_data);
}

} // end namespace profugus

#endif // utils_Metaclass_t_hh

//---------------------------------------------------------------------------//
//                 end of Metaclass.t.hh
//---------------------------------------------------------------------------//
