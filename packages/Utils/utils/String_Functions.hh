//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/String_Functions.hh
 * \author Seth R Johnson
 * \date   2012/01/31
 * \brief  String-related functions.
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_String_Functions_hh
#define utils_String_Functions_hh

#include <string>
#include <sstream>

#include "harness/DBC.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \fn to_lower
 * \brief Returns a lowercased std::string given an ASCII std::string input
 *
 */
/*!
 * \example utils/test/tstString_Functions.cc
 *
 * String functions test
 */
//===========================================================================//

std::string to_lower(const std::string& orig_string);

//===========================================================================//
/*!
 * \fn join
 * \brief similar to python's "str.join" method
 *
 * This joins all given elements, inserting conjunction betwen them.
 */
//===========================================================================//

template<class InputIterator>
std::string join(
        InputIterator first,
        InputIterator last,
        const std::string& conjunction)
{
    std::ostringstream result;
    InputIterator it = first;

    // First element is not preceded by a conjunction
    result << *it++;
    // Join the rest
    while (it != last)
        result << conjunction << *it++;

    return result.str();
}

} // end namespace profugus

#endif // utils_String_Functions_hh

//---------------------------------------------------------------------------//
//              end of utils/String_Functions.hh
//---------------------------------------------------------------------------//
