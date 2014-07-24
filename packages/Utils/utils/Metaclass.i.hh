//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Metaclass.i.hh
 * \author Seth R Johnson
 * \date   Fri Aug 30 20:59:34 2013
 * \brief  Member definitions of class Metaclass.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Metaclass_i_hh
#define utils_Metaclass_i_hh

namespace profugus
{
//---------------------------------------------------------------------------//
// STATIC METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Add a new member to the class, returning its entry ID.
 *
 * Note that the name is solely for user reference; because of the performance
 * penalty of doing string lookups, we never search over strings.
 *
 * This takes into account alignment of items in the structure. For example, if
 * you add a 1-byte char and then a 4-byte float, the "packed" alignment may
 * look like <pre>
       c f f f f
 * </pre>
 * for a total of 5 bytes, but reading the float out of that buffer repeatedly
 * is on many architectures more expensive than aligning it to a 4-byte
 * boundary: <pre>
       c _ _ _ f f f f
 * </pre>
 */
template<typename I>
template<typename T>
unsigned int Metaclass<I>::new_member(
    const std::string& name,
    SP_Member_Manager  manager)
{
    Require(name.size() > 0);
    Require(manager);
    Remember(unsigned int orig_storage_size = d_storage_size);
    Remember(unsigned int orig_size = d_member_md.size());

#if UTILS_DBC > 0
    Insist(Metaclass::d_num_instances == 0,
           "new_member can only be called when no instances exist.");
#endif

    // >>> Add a new entry
    d_member_md.resize(d_member_md.size() + 1);

    Member_Metadata& md = d_member_md.back();

    md.name = name;
#ifdef REQUIRE_ON
    md.type_string = typeid(T).name();
#endif
    md.member_manager = manager;

    // >>> Determine the offset of this data item

    unsigned int aligned_size = sizeof(T);
#if ((__cplusplus >= 201103L) && (!defined(__INTEL_COMPILER)))
    aligned_size = alignof(T);
#endif

    md.data_offset = ((d_storage_size + aligned_size - 1) / aligned_size)
                     * aligned_size;
    Check(md.data_offset % aligned_size == 0);
    Check(md.data_offset >= d_storage_size);

    // Now update the accumulated size of this class
    d_storage_size = md.data_offset + sizeof(T);

    Ensure(d_storage_size > orig_storage_size);
    Ensure(d_member_md.size() > orig_size);
    return d_member_md.size() - 1;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access data for a given entry ID
 */
template<typename I>
template<typename T>
inline T& Metaclass<I>::access(unsigned int member_id)
{
    Require(d_data);
    Require(is_valid_access<T>(member_id));

    // Pointer to the start of the requested entry
    char* data = d_data + d_member_md[member_id].data_offset;

    // Make sure the back end of data doesn't point past the end of allocated
    // memory
    Ensure(data + sizeof(T) <= data + storage_size());
    return reinterpret_cast<T&>(*data);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access data for a given entry ID (const)
 */
template<typename I>
template<typename T>
inline const T& Metaclass<I>::access(unsigned int member_id) const
{
    Require(d_data);
    Require(is_valid_access<T>(member_id));

    // Pointer to the start of the requested entry
    const char* data = d_data + d_member_md[member_id].data_offset;

    // Make sure the back end of data doesn't point past the end of allocated
    // memory
    Ensure(data + sizeof(T) <= data + storage_size());
    return reinterpret_cast<const T&>(*data);
}

//---------------------------------------------------------------------------//

template<typename I>
template<typename T>
inline bool Metaclass<I>::is_valid_access(unsigned int member_id) const
{
    if (member_id >= size())
        return false;

    bool type_equal = true;
#ifdef REQUIRE_ON
    type_equal = (typeid(T).name() == d_member_md[member_id].type_string);
    Validate(type_equal,
             "Tried to access member '" << d_member_md[member_id].name
             << "' (index " << member_id << ") with incorrect type.");
#endif
    return type_equal;
}

} // end namespace profugus

#endif // utils_Metaclass_i_hh

//---------------------------------------------------------------------------//
//                 end of Metaclass.i.hh
//---------------------------------------------------------------------------//
