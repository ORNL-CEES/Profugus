//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/Demangle_Type.cc
 * \author Seth R Johnson
 * \date   Wed Oct 16 08:23:54 2013
 * \brief  Demangle_Type member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Demangle_Type.hh"

#ifdef __GNUG__
#include <cstdlib>
#include <memory>
#include <cxxabi.h>
#endif // __GNUG__

namespace nemesis
{

std::string Demangle_Type::demangle(const char* typeid_name)
{
#ifdef __GNUG__
    int status = -1;
    // Return a null-terminated string allocated with malloc
    char* demangled = abi::__cxa_demangle(typeid_name, NULL, NULL, &status);

    // Copy the C string to a STL string if successful, or the mangled name if
    // not
    std::string result(status == 0 ? demangled : typeid_name);

    // Free the returned memory
    free(demangled);
#else // __GNUG__
    std::string result(typeid_name);
#endif // __GNUG__

    return result;
}

} // end namespace nemesis

//---------------------------------------------------------------------------//
//                 end of Demangle_Type.cc
//---------------------------------------------------------------------------//
