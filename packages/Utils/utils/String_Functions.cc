//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/String_Functions.cc
 * \author Seth R. Johnson
 * \date   Sat Sep 03 18:55:17 2011
 * \brief  Member definitions of String_Functions classes.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <string>
#include "String_Functions.hh"

namespace profugus
{

std::string to_lower(const std::string& orig_string)
{
    std::string outString(orig_string);
    std::transform(outString.begin(), outString.end(),
            outString.begin(), ::tolower);

    return outString;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of String_Functions.cc
//---------------------------------------------------------------------------//
