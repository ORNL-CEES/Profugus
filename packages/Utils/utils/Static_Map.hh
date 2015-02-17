//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Static_Map.hh
 * \author Seth R Johnson
 * \date   Tue Dec 17 21:11:12 2013
 * \brief  Static_Map class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Static_Map_hh
#define utils_Static_Map_hh

#include <algorithm>
#include <vector>
#include <utility>

#include "Teuchos_RCP.hpp"
#include "harness/DBC.hh"

#include "Default_Hash.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Static_Map
 * \brief Defines a fixed-size hash table
 *
 * This class defines a hash table in which the size is \b fixed after
 * construction by a call to \c complete.  It does not support dynamic size;
 * insertion and deletion are not allowed after the table is
 * "completed". Internally, the table may store more buckets than objects
 * depending on the requirements of the defined Hash_Function.
 *
 * This hash table supports iterating over its contents after completion.
 * Items may be iterated upon while building the hash table, but I don't know
 * why you'd want to do this. It also invalidates the iterators.
 *
 * The hash table performs conflict resolution using \e seperate-chaining.
 *
 * The hash table can be built with different \c Hash_Functions that are
 * included as a template parameter; it defaults to a simple modulo hash
 * function.
 *
 * The default hash function (modulo) has the fewest collisions if the numbers
 * are all consecutive. If this requirement is fulfilled, accessing should be
 * nearly as efficient as pulling from an array.
 *
 * \sa Int_Mod_Hash_Function
 */
/*!
 * \example utils/test/tstStatic_Map.cc
 *
 * Test of Static_Map.
 */
//===========================================================================//

template<class Key, class T, class Hash=profugus::Default_Hash<Key> >
class Static_Map
{
    typedef Static_Map<Key, T, Hash> This;

  public:
    //@{
    //! Template typedefs.
    typedef Key  key_type;
    typedef T    mapped_type;
    typedef Hash hasher;
    //@}

    //! Results returned by the map (key/value pair)
    typedef std::pair<Key, T> value_type;
    //! Vector of keys and values (implementation detail)
    typedef std::vector<value_type> Container_t;

    //@{
    //! More container typedefs.
    typedef typename Container_t::iterator        iterator;
    typedef typename Container_t::const_iterator  const_iterator;
    typedef typename Container_t::size_type       size_type;
    typedef typename Container_t::difference_type difference_type;
    typedef typename Container_t::reference       reference;
    typedef typename Container_t::const_reference const_reference;
    //@}

  private:
    // >>> DATA and TYPES

    // Key-value item stored in the list.
    typedef value_type Item;
    // Begin/end indices
    typedef std::pair<size_type, size_type> Pair_Size;
    // Vector of begin/end indices
    typedef std::vector<Pair_Size> Vec_Pair_Size;
    // Hash function
    typedef Teuchos::RCP<hasher> RCP_Hash;

    //! Hash function.
    RCP_Hash d_hash;

    // Vector of actual keys and values: before "complete" is called, unsorted;
    // afterward, sorted according to hash of key. "Buckets" are contiguous
    // subsets of this list.
    Container_t d_items;

    // Buckets: begin/end incices
    Vec_Pair_Size d_buckets;

  public:
    // >>> CONSTRUCTORS

    //! Default constructor
    Static_Map()
    {
        this->clear();
    }

    //! Construct from key/value iterators, e.g. a map
    template<typename InputIterator>
    Static_Map(InputIterator first, InputIterator last)
    {
        this->clear();
        this->insert(first, last);
    }

    //! Copy constructor
    Static_Map(const This& rhs)
    {
        *this = rhs;
        ENSURE(completed() == rhs.completed());
        ENSURE(size() == rhs.size());
    }

    // Assignment operator
    This& operator= (const This& rhs)
    {
        d_hash    = rhs.d_hash;
        d_items   = rhs.d_items;
        d_buckets = rhs.d_buckets;

        ENSURE(completed() == rhs.completed());
        ENSURE(size() == rhs.size());
        return *this;
    }

    // >>> CONSTRUCTION

    // Insert a unique item.
    inline void insert(const value_type& value);

    // Insert a bunch of items as key/value pairs
    template<typename InputIterator>
    void insert(InputIterator iter, InputIterator last)
    {
        for (; iter != last; ++iter)
        {
            insert(*iter);
        }
    }

    // Construct the hash table.
    void complete();

    // Get the value (after complete has been called).
    inline T& operator[](key_type key);

    // Get the value (after complete has been called).
    inline const T& operator[](key_type key) const;

    //@{
    //! Get the value (after it has been inserted).
    inline T& at(key_type key) { return operator[](key); }
    inline const T& at(key_type key) const { return operator[](key); }
    //@}

    // Find an element in the map
    const_iterator find(key_type key) const;

    /*! \brief Return the number of items matching 'key'
     *
     * This will always return 0 (not present in table) or 1 (present); mostly
     * for STL compatibility.
     */
    size_type count(key_type key) const
    {
        REQUIRE(this->completed());
        return (exists(key) ? 1u : 0u);
    }

    // Search to see if an object is in the hash table.
    inline bool exists(key_type key) const;

    // Remove all entries from the hash table and set to "incomplete"
    void clear()
    {
        d_hash = RCP_Hash();
        d_items.clear();
        d_buckets.clear();

        ENSURE(!completed());
        ENSURE(empty());
        ENSURE(d_hash.is_null());
    }

    // Remove all entries from the hash table and set to "incomplete"
    void swap(This& rhs)
    {
        std::swap(d_hash   , rhs.d_hash);
        std::swap(d_items  , rhs.d_items);
        std::swap(d_buckets, rhs.d_buckets);
    }

    // >>> ACCESSORS

    //! Whether we've finished being filled and are ready to access
    bool completed() const
    {
        return !d_buckets.empty();
    }

    //! Whether we have no objects
    bool empty() const { return d_items.empty(); }

    //! Number of inserted objects.
    size_type size() const { return d_items.size(); }

    // >>> ITERATORS

    //@{
    //! Iterator access to items before or after completion
    const_iterator begin() const { return d_items.begin(); }
    const_iterator end()   const { return d_items.end(); }
    //@}

    // >>> BUCKET INTERFACE

    //! Number of buckets required to store the objects.
    size_type bucket_count() const
    {
        REQUIRE(completed());
        return d_buckets.size();
    }

    //! Number of items in a bucket.
    inline size_type bucket_size(size_type bucket_id) const;

    // Return the maximum number of items in any bucket.
    inline size_type max_items_per_bucket() const;

    //! Return the average number of elements per bucket
    float load_factor() const
    {
        return static_cast<float>(size()) / bucket_count();
    }

    // >>> HASH INTERFACE
    const hasher& hash_function() const
    {
        REQUIRE(completed());
        ENSURE(!d_hash.is_null());
        return *d_hash;
    }

  private:
    // >>> IMPLEMENTATION

    // Sort first on hash value, then key value
    struct Hash_Comparator
    {
        const Hash& d_hash;

        Hash_Comparator(const Hash& hash) : d_hash(hash) { /* * */ }

        bool operator() (const value_type& left, const value_type& right) const
        {
            const size_type left_hash  = d_hash.hash(left.first);
            const size_type right_hash = d_hash.hash(right.first);

            if (left_hash == right_hash)
            {
                // sort by key inside the hash
                return left.first < right.first;
            }

            return left_hash < right_hash;
        }
    };
};

template<class K, class T, class H>
inline void swap(Static_Map<K,T,H>& a, Static_Map<K,T,H>& b)
{
    a.swap(b);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
// TEMPLATE AND INLINE FUNCTION DEFS
//---------------------------------------------------------------------------//

#include "Static_Map.i.hh"

#endif // utils_Static_Map_hh

//---------------------------------------------------------------------------//
//              end of utils/Static_Map.hh
//---------------------------------------------------------------------------//
