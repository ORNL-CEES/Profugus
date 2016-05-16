//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/Member_Manager.hh
 * \author Seth R Johnson
 * \date   Wed Oct 16 07:36:23 2013
 * \brief  Member_Manager and subclass class definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Member_Manager_hh
#define Utils_utils_Member_Manager_hh

#include <cstring>
#include <cstdlib>
#include <vector>

#include "harness/DBC.hh"

namespace profugus
{

namespace metaclass
{

//===========================================================================//
/*!
 * \class Member_Manager
 * \brief Manage creation/destruction of a Metaclass member
 *
 * This virtual base class must have all its methods implemented by whatever
 * subclasses it.
 *
 * \sa Metaclass
 */
//===========================================================================//

class Member_Manager
{
  public:
    typedef std::size_t size_type;

  public:
    virtual ~Member_Manager();

    /*!
     * \brief Interface for constructor
     *
     * This is given the already-allocated memory (using malloc), so it has to
     * do a "placement new" construction if the object has a nontrivial
     * constructor. Otherwise, it can use memset.
     *
     * \param[in,out] data Address of memory to initialize
     */
    virtual void construct(void* data) = 0;

    /*!
     * \brief Interface for copy constructor
     *
     * This is given the already-allocated memory (using malloc), so it has to
     * do a "placement new" construction if the object has a nontrivial
     * constructor.
     *
     * \param[in,out] data Address of memory to initialize
     * \param[in]     rhs_data Address of object to copy
     */
    virtual void copy_construct(void* data, const void* rhs_data) = 0;

    /*!
     * \brief Construct by unpacking from a buffer
     *
     * \param[in] data Address of object to initialize
     * \param[in] buffer Buffer from which to unpack
     *
     * \return Pointer to the "end" of the data just unpacked
     */
    virtual const char* unpack_construct(void* data, const char* buffer) = 0;

    /*!
     * \brief Interface for assignment
     *
     * \param[in,out] data     Address of initialized object to assign to
     * \param[in]     rhs_data Address of object to copy
     */
    virtual void assign(void* data, const void* rhs_data) = 0;

    /*!
     * \brief Interface for destructor
     *
     * This is meant to call the destructor, not deallocate.
     *
     * \param[in,out] data     Address of initialized object to destroy
     */
    virtual void destroy(void* data) = 0;

    /*!
     * \brief Pack the data
     *
     * \param[in] data Address of object to pack
     * \param[in] buffer Buffer into which to pack
     *
     * \return Pointer to the "end" of the data just packed in
     */
    virtual char* pack(const void* data, char* buffer) = 0;

    //! Packed size
    virtual size_type packed_size(const void* data) const = 0;
};

//===========================================================================//
/*!
 * \brief Manage construction/deletion PODs
 */
//===========================================================================//

template<class T>
class Member_Manager_POD : public Member_Manager
{
  public:
    typedef T Data_t;

  public:
    //! If debugging, clear the memory.
    void construct(void* data)
    {
#ifdef CHECK_ON
        std::memset(data, 0, sizeof(Data_t));
#endif
    }

    //! Use the copy constructor to initialize already allocated memory
    void copy_construct(void* data, const void* rhs_data)
    {
        std::memcpy(data, rhs_data, sizeof(Data_t));
    }

    //! Create from a packed buffer
    const char* unpack_construct(void* data, const char* buffer)
    {
        std::memcpy(data, buffer, sizeof(Data_t));
        return buffer + sizeof(Data_t);
    }

    //! Use the copy constructor to initialize already allocated memory
    void assign(void* data, const void* rhs_data)
    {
        std::memcpy(data, rhs_data, sizeof(Data_t));
    }

    //! Destruction needs no special treatment
    void destroy(void* ) { /* * */ }

    //! Copy data into the pack buffer
    char* pack(const void* data, char* buffer)
    {
        std::memcpy(buffer, data, sizeof(Data_t));
        return buffer + sizeof(Data_t);
    }

    //! Size required for packing is just size of this type
    virtual size_type packed_size(const void*) const { return sizeof(Data_t); }
};

//===========================================================================//
/*!
 * \brief Manage construction/deletion of STL vectors
 */
//===========================================================================//

template<class T>
class Member_Manager_Vector : public Member_Manager
{
  public:
    typedef T                       value_type;
    typedef std::vector<value_type> Vector_t;

  public:

    //! Call the constructor on already allocated memory
    void construct(void* data)
    {
        new (data) Vector_t();
    }

    //! Use the copy constructor to initialize already allocated memory
    virtual void copy_construct(void* data, const void* rhs_data)
    {
        const Vector_t* rhs = reinterpret_cast<const Vector_t*>(rhs_data);
        // Calling this constructor then tells the STL vector to allocate space
        // for all the entries, and to copy those.
        new (data) Vector_t(*rhs);
    }

    virtual const char* unpack_construct(void* data, const char* buffer)
    {
        size_type num_items;
        std::memcpy(&num_items, buffer, sizeof(size_type));
        buffer += sizeof(size_type);

        // Call STL vector constructor with iterators
        const value_type* begin = reinterpret_cast<const value_type*>(buffer);
        new (data) Vector_t(begin, begin + num_items);
        buffer += num_items * sizeof(value_type);

        return buffer;
    }

    //! Use operator= for assignment
    virtual void assign(void* lhs_data, const void* rhs_data)
    {
        Vector_t* lhs       = reinterpret_cast<Vector_t*>(lhs_data);
        const Vector_t* rhs = reinterpret_cast<const Vector_t*>(rhs_data);

        *lhs = *rhs;
    }

    //! Call the destructor
    void destroy(void* data)
    {
        Vector_t* vec = reinterpret_cast<Vector_t*>(data);

        // Call destructor
#ifdef __clang__
        vec->~vector();
#else
        // Some older compilers may still require this?
        vec->template ~Vector_t();
#endif
    }

    //! Pack all items from the vector into the buffer
    virtual char* pack(const void* data, char* buffer)
    {
        const Vector_t& vec = *reinterpret_cast<const Vector_t*>(data);

        // Copy the number of elements
        const size_type num_items = vec.size();
        std::memcpy(buffer, &num_items, sizeof(size_type));
        buffer += sizeof(size_type);

        // Copy the elements
        std::memcpy(buffer, &vec[0], num_items * sizeof(value_type));
        buffer += num_items * sizeof(value_type);

        return buffer;
    }

    virtual size_type packed_size(const void* data) const
    {
        const Vector_t& vec = *reinterpret_cast<const Vector_t*>(data);
        return sizeof(size_type) + vec.size() * sizeof(value_type);
    }
};

} // end namespace profugus::metaclass

} // end namespace profugus

#endif // Utils_utils_Member_Manager_hh

//---------------------------------------------------------------------------//
//                 end of Member_Manager.hh
//---------------------------------------------------------------------------//
