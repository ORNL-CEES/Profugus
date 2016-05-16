//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/String_Functions.hh
 * \author Seth R Johnson
 * \date   2012/01/31
 * \brief  String-related functions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_String_Functions_hh
#define Utils_utils_String_Functions_hh

#include <string>
#include <sstream>

#include "harness/DBC.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \fn lower
 * \brief Returns a lowercased std::string given an ASCII std::string input
 *
 */
/*!
 * \example utils/test/tstString_Functions.cc
 *
 * String functions test
 */
//===========================================================================//

std::string lower(const std::string& orig_string);

//===========================================================================//
/*!
 * \fn startswith
 * \brief Whether the string starts with another string.
 */
//===========================================================================//

bool startswith(const std::string& main_string, const std::string& suffix);

//===========================================================================//
/*!
 * \fn endswith
 * \brief Whether the string ends with another string.
 */
//===========================================================================//

bool endswith(const std::string& main_string, const std::string& suffix);

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
    if (it != last)
        result << *it++;

    // Join the rest
    while (it != last)
        result << conjunction << *it++;

    return result.str();
}

//===========================================================================//
/*!
 * \fn rstrip
 * \brief Removes trailing whitespace from a string.
 */
//===========================================================================//

std::string rstrip(const std::string &s);

//===========================================================================//
/*!
 * \fn max_length
 * \brief Find the maximum length of a series of strings.
 *
 * This is useful for creating tables.
 */
//===========================================================================//
template<class InputIterator>
std::size_t max_length(
        InputIterator first,
        InputIterator last)
{
    size_t result = 0;
    while (first != last)
    {
        auto current = first->size();
        if (current > result)
            result = current;
        ++first;
    }
    return result;
}

} // end namespace profugus

#endif // Utils_utils_String_Functions_hh

//---------------------------------------------------------------------------//
//              end of utils/String_Functions.hh
//---------------------------------------------------------------------------//
