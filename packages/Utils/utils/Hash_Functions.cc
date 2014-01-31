//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Hash_Functions.cc
 * \author Thomas M. Evans
 * \date   Sat Sep 03 18:55:17 2011
 * \brief  Member definitions of Hash_Functions classes.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Hash_Functions.hh"

//===========================================================================//
// PRIME NUMBER TABLE
//===========================================================================//

namespace
{

//! Table of prime numbers.
def::size_type primes[] = {251,         // 0
                           509,         // 1
                           1021,        // 2
                           2039,        // 3
                           4093,        // 4
                           8191,        // 5
                           16381,       // 6
                           32749,       // 7
                           65521,       // 8
                           131071,      // 9
                           262139,      // 10
                           524287,      // 11
                           1048573,     // 12
                           2097143,     // 13
                           4194301,     // 14
                           8388593,     // 15
                           16777213,    // 16
                           33554393,    // 17
                           67108859,    // 18
                           134217689,   // 19
                           268435399,   // 20
                           536870909,   // 21
                           1073741789,  // 22
                           2147483647}; // 23

//! Number of primes in the table.
const def::size_type num_primes = 24;

}

//===========================================================================//
// HASH_FUNCTIONS DEFINITIONS
//===========================================================================//

namespace profugus
{

//---------------------------------------------------------------------------//
// INT_MOD_HASH_FUNCTION DEFINITIONS
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Int_Mod_Hash_Function::Int_Mod_Hash_Function(size_type num_objects)
{
    // find the prime number > num_objects using the prime number table

    size_type *M = std::lower_bound(&primes[0],
                                    &primes[num_primes], num_objects);
    Check (M != &primes[num_primes]);

    // assign M
    d_M = *M;

    Ensure (d_M >= num_objects);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Hash_Functions.cc
//---------------------------------------------------------------------------//
