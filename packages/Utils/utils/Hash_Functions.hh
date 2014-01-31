//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Hash_Functions.hh
 * \author Thomas M. Evans
 * \date   Sat Sep 03 18:52:06 2011
 * \brief  Hash_Functions class definition.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Hash_Functions_hh
#define utils_Hash_Functions_hh

#include <algorithm>

#include "harness/DBC.hh"
#include "Definitions.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Int_Mod_Hash_Function
 * \brief Defines the modular hash function \f$k\bmod M\f$ for unsigned integer
 * keys.
 *
 * This class defines the modular hash function
 * \f[
   h = k \bmod M
 * \f]
 * where \e k is the integer key and \e M is a prime number \c >= the number
 * of objects stored in the hash table.
 *
 * \sa Static_Hash_Table
 */
/*!
 * \example utils/test/tstHash_Functions.cc
 *
 * Hash functions test
 */
//===========================================================================//

class Int_Mod_Hash_Function
{
  public:
    //@{
    //! Required typedefs.
    typedef def::size_type size_type;
    typedef def::size_type key_type;
    //@}

  private:
    // >>> DATA

    // Number of buckets this function requires to store the requested number
    // of objects (a prime > number of objects).
    size_type d_M;

  public:
    // Constructor.
    explicit Int_Mod_Hash_Function(size_type num_objects);

    // Number of buckets required to store the number of objects in a table.
    size_type num_buckets() const { return d_M; }

    // Hash function.
    size_type hash(key_type k) const { return k % d_M; }
};

} // end namespace profugus

#endif // utils_Hash_Functions_hh

//---------------------------------------------------------------------------//
//              end of utils/Hash_Functions.hh
//---------------------------------------------------------------------------//
