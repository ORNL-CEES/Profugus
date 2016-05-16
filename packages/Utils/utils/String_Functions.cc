//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/String_Functions.cc
 * \author Seth R. Johnson
 * \date   Sat Sep 03 18:55:17 2011
 * \brief  Member definitions of String_Functions classes.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//


#include "String_Functions.hh"

#include <algorithm>
#include <string>
#include <cctype>


namespace profugus
{
//---------------------------------------------------------------------------//
std::string lower(const std::string& orig_string)
{
    std::string outString(orig_string);
    std::transform(outString.begin(), outString.end(),
            outString.begin(), ::tolower);

    return outString;
}

//---------------------------------------------------------------------------//
std::string rstrip(const std::string &s)
{
    // Copy the string to get a working function
    std::string ws = s;

    // Erase the trailing whitespace
    ws.erase(std::find_if(ws.rbegin(), ws.rend(),
                          [](char c) { return !std::isspace(c); }).base(),
             ws.end());
    return ws;
}

//---------------------------------------------------------------------------//
bool startswith(const std::string& main_string, const std::string& prefix)
{
    if (main_string.size() < prefix.size())
        return false;

    return std::equal(main_string.begin(),
                      main_string.begin() + prefix.size(),
                      prefix.begin());
}

//---------------------------------------------------------------------------//
bool endswith(const std::string& main_string, const std::string& suffix)
{
    if (main_string.size() < suffix.size())
        return false;

    return std::equal(main_string.end() - suffix.size(), main_string.end(),
                      suffix.begin());
}


} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of String_Functions.cc
//---------------------------------------------------------------------------//
