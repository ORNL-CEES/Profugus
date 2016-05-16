//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/Default_Hash.hh
 * \author Seth R Johnson
 * \date   Fri Dec 06 10:06:08 2013
 * \brief  Default_Hash class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Default_Hash_hh
#define Utils_utils_Default_Hash_hh

#include "Definitions.hh"
#include "Hash_Functions.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \struct Default_Hash
 * \brief Chooses a default hash function based on key type.
 *
 * The structure that gets created must have a
 */
//===========================================================================//

template<class Key>
struct Default_Hash;

//===========================================================================//
template<>
struct Default_Hash<def::size_type> : public Int_Mod_Hash_Function
{
    explicit Default_Hash(size_type num_objects)
        : Int_Mod_Hash_Function(num_objects)
    { /* * */ }
};

//===========================================================================//
template<>
struct Default_Hash<int> : public Int_Mod_Hash_Function
{
    explicit Default_Hash(size_type num_objects)
        : Int_Mod_Hash_Function(num_objects)
    { /* * */ }
};

//===========================================================================//
template<>
struct Default_Hash<unsigned int> : public Int_Mod_Hash_Function
{
    explicit Default_Hash(size_type num_objects)
        : Int_Mod_Hash_Function(num_objects)
    { /* * */ }
};

//===========================================================================//

} // end namespace profugus

#endif // Utils_utils_Default_Hash_hh

//---------------------------------------------------------------------------//
//                 end of Default_Hash.hh
//---------------------------------------------------------------------------//
