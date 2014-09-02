//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Static_Map.i.hh
 * \author Seth R Johnson
 * \date   Tue Dec 17 21:11:08 2013
 * \brief  Member definitions of class Static_Map.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Static_Map_i_hh
#define utils_Static_Map_i_hh

namespace profugus
{

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Insert an object into the map
 *
 * \pre complete() has not been called
 *
 * \warning This operation will not fail immediately if you insert a duplicate
 * key; it will only fail when complete() is called.
 */
template<class Key, class T, class Hash>
void Static_Map<Key,T,Hash>::insert(const value_type& value)
{
    REQUIRE(!completed());

    d_items.push_back(value);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Change hash table from "building" phase to accessible phase
 *
 * This step is an O(log(N)) sort plus an O(N) verification/construction.
 *
 * \note This routine sorts the item list according to hash key, then loops
 * through the list again to determine which indices in the list correspond to
 * what "bucket" in the hash table.
 * The result is that d_buckets[i] corresponds to the begin/end indices inside
 * d_items for hash value i.
 */
template<class Key, class T, class Hash>
void Static_Map<Key,T,Hash>::complete()
{
    REQUIRE(!completed());

    // Create the hash using the number of objects in the hash table
    d_hash = Teuchos::rcp(new Hash(d_items.size()));

    // Create empty buckets (begin == end)
    CHECK(!d_hash.is_null());
    d_buckets.assign(d_hash->num_buckets(), Pair_Size(0, 0));

    if (d_items.empty())
    {
        ENSURE(completed());
        ENSURE(empty());
        return;
    }

    // Sort the key/value pairs according to hash
    std::sort(d_items.begin(), d_items.end(), Hash_Comparator(*d_hash));

    // Put them into consecutive buckets according to hash, while checking for
    // uniqueness by comparing against previous key

    key_type  prev_key   = d_items[0].first;
    size_type prev_hash  = d_hash->hash(prev_key);

    CHECK(prev_hash < d_buckets.size());

    for (size_type i = 1; i < d_items.size(); ++i)
    {
        const key_type this_key = d_items[i].first;
        VALIDATE(this_key != prev_key, "Duplicate entry " << this_key
                << " in static hash table.");

        const size_type this_hash = d_hash->hash(this_key);
        if (this_hash == prev_hash)
            continue;

        CHECK(this_hash > prev_hash);

        // If different hash, the last bucket has ended...
        d_buckets[prev_hash].second = i;
        // ...and the new bucket has begun
        d_buckets[this_hash].first  = i;

        prev_hash = this_hash;
    }
    // Add the final bucket 'end' index
    d_buckets[prev_hash].second = d_items.size();

    ENSURE(completed());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the value for a key.
 *
 * The \c (key,value) pair must have been inserted in the database.  For
 * efficiency, we do not search to see if the key is not already defined.
 *
 * \pre \c (key,value) pair is already inserted in the database
 */
template<class Key, class T, class Hash>
const T& Static_Map<Key,T,Hash>::operator[](key_type key) const
{
    REQUIRE(completed());
    REQUIRE(exists(key));

    const_iterator iter = d_items.begin() + d_buckets[d_hash->hash(key)].first;
    do
    {
        if (key == iter->first)
            return iter->second;

        ++iter;
    } while (true);

    // The following line should never be called and will be compiled out
    // (at least as tested on GCC 4.7
    return iter->second;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the mutable value for a key.
 *
 * The \c (key,value) pair must have been inserted in the database.  For
 * efficiency, we do not search to see if the key is not already defined.
 *
 * \pre \c complete()  has already been called
 * \pre \c (key,value) pair is already inserted in the database
 */
template<class Key, class T, class Hash>
T& Static_Map<Key,T,Hash>::operator[](key_type key)
{
    REQUIRE(completed());
    REQUIRE(exists(key));

    iterator iter = d_items.begin() + d_buckets[d_hash->hash(key)].first;
    do
    {
        if (key == iter->first)
            return iter->second;

        ++iter;
    } while (true);

    // The following line should never be called and will be compiled out
    // (at least as tested on GCC 4.7
    return iter->second;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find an item in the map, or return end() if not found
 *
 * This allows you to quickly check for both existence or access it, with a
 * single hash evaluation and traversal. It also provides STL compliance.
 *
 * \pre \c complete()  has already been called
 */
template<class Key, class T, class Hash>
typename Static_Map<Key,T,Hash>::const_iterator
Static_Map<Key,T,Hash>::find(key_type key) const
{
    REQUIRE(completed());

    const Pair_Size begend = d_buckets[d_hash->hash(key)];

    const_iterator iter = d_items.begin() + begend.first;
    const_iterator last = d_items.begin() + begend.second;

    for (; iter != last; ++iter)
    {
        if (key == iter->first)
            return iter;
    }

    // Not found
    return end();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Search the hash table to see if a key has been inserted.
 *
 * \param key key to search the database for
 * \return true if key has been inserted, false otherwise
 */
template<class Key, class T, class Hash>
bool Static_Map<Key,T,Hash>::exists(key_type key) const
{
    return find(key) != end();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the number of items in a bucket.
 */
template<class Key, class T, class Hash>
typename Static_Map<Key,T,Hash>::size_type
Static_Map<Key,T,Hash>::bucket_size(size_type bucket_id) const
{
    REQUIRE(completed());
    REQUIRE(bucket_id < d_buckets.size());

    const Pair_Size& begend = d_buckets[bucket_id];

    return begend.second - begend.first;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the maximum number of items in any one bucket.
 */
template<class Key, class T, class Hash>
typename Static_Map<Key,T,Hash>::size_type
Static_Map<Key,T,Hash>::max_items_per_bucket() const
{
    // max items
    size_type n = 0;

    // loop through buckets and calculate the maximum number of items in any
    // one bucket
    for (typename Vec_Pair_Size::const_iterator iter = d_buckets.begin();
         iter != d_buckets.end(); ++iter)
    {
        CHECK(iter->second >= iter->first);
        n = std::max(n, iter->second - iter->first);
    }
    return n;
}

} // end namespace profugus

#endif // utils_Static_Map_i_hh

//---------------------------------------------------------------------------//
//              end of utils/Static_Map.i.hh
//---------------------------------------------------------------------------//
