//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Metaclass.hh
 * \author Seth R Johnson
 * \date   Fri Aug 30 20:59:34 2013
 * \brief  Metaclass class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Metaclass_hh
#define utils_Metaclass_hh

#include <vector>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <memory>

#include "harness/DBC.hh"
#include "Member_Manager.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Metaclass
 * \brief Manage and instantiate runtime "classes".
 *
 * This is used for creating an efficient struct accessor at run-time.
 *
 * \tparam Instantiator Class used to make distinct compile-time structs.
 */
/*!
 * \example utils/test/tstMetaclass.cc
 *
 * Test of Metaclass.
 */
//===========================================================================//

template<class Instantiator>
class Metaclass
{
    typedef Metaclass<Instantiator> This;

  public:
    //@{
    //! Typedefs.
    typedef std::shared_ptr<metaclass::Member_Manager> SP_Member_Manager;
    typedef std::size_t                                size_type;
    //@}

  private:
    // >>> Underlying data storage (size is known and static)
    char* d_data;

  public:
    // >>> CONSTRUCTION

    // Construct
    Metaclass();

    // Copy construct
    Metaclass(const This& rhs);

    // Construct from packed data
    Metaclass(const char* const buffer, unsigned int size);

    // Clean up on destruction
    ~Metaclass();

    // Assignment operator
    This& operator=(const This& rhs);

    // >>> PACKING

    // Get the pickled size of this struct in bytes (nontrivial)
    size_type packed_size() const;

    // Pack into a properly-sized buffer
    char* pack(char* buffer) const;

    // Access a given entry in the dynamic struct.
    template<typename T>
    inline T& access(unsigned int member_id);

    // Access a given entry in the dynamic struct. (const)
    template<typename T>
    inline const T& access(unsigned int member_id) const;

  private:
    // >>> IMPLEMENTATION

    // Validate type, only used if REQUIRE_ON
    template<typename T>
    inline bool is_valid_access(unsigned int member_id) const;

    // Allocate the memory for an instance
    void alloc();

  public:
    // >>> STATIC METHODS

    // Add a new member to the class, returning its entry ID.
    template<typename T>
    static unsigned int new_member(const std::string& name,
                                   SP_Member_Manager  manager);

    /*!
     * \brief Add a new plain-old-data member to the class
     *
     * This is a convenience function for adding ints, simple structs,
     * Vector_Lite, etc.
     *
     * \return the new entry ID.
     */
    template<typename T>
    static unsigned int new_pod_member(const std::string& name)
    {
        SP_Member_Manager manager(
            std::make_shared<metaclass::Member_Manager_POD<T>>());

        return new_member<T>(name, manager);
    }

    /*!
     * \brief Add a new vector member to the class
     *
     * This is a convenience function for adding simple vectors only.
     *
     * \tparam Vector_T the fully qualified vector type
     * \return the new entry ID.
     */
    template<typename Vector_T>
    static unsigned int new_vec_member(const std::string& name)
    {
        typedef typename Vector_T::value_type value_type;

        SP_Member_Manager manager(
            std::make_shared<metaclass::Member_Manager_Vector<value_type>>());

        return new_member<Vector_T>(name, manager);
    }

    // Access the name of an entry in the dynamic struct
    static std::string name(unsigned int member_id);

    // Reset all entries
    static void reset();

    //! Number of entries
    static size_type size() { return d_member_md.size(); }

    //! Amount of data taken up by each instance
    static size_type storage_size() { return d_storage_size; }

  private:
    //! This data struct is for the static array that determines entry offsets
    struct Member_Metadata
    {
        // Descriptive string
        std::string name;

#ifdef REQUIRE_ON
        // Class name from typeid(T).name()
        std::string type_string;
#endif

        // Memory management
        SP_Member_Manager member_manager;

        // Index into the data entry
        unsigned int data_offset;
    };

    // Offset info for each extra data type
    typedef std::vector<Member_Metadata> Vec_Member_Md;

  private:
    // >>> STATIC DATA

    // Metadata about members
    static Vec_Member_Md d_member_md;

    // Accumulated member data size (bytes)
    static unsigned int d_storage_size;

    // Total number of extant instances, for error checking
    static unsigned int d_num_instances;
};

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//

#include "Metaclass.i.hh"

//---------------------------------------------------------------------------//

#endif // utils_Metaclass_hh

//---------------------------------------------------------------------------//
//                 end of Metaclass.hh
//---------------------------------------------------------------------------//
