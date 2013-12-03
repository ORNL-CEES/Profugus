//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/Release.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 11:50:34 2008
 * \brief  Release function for the harness library
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef harness_Release_hh
#define harness_Release_hh

#include <string>

//===========================================================================//
/*!
 * \namespace nemesis
 *
 * \brief Namespace that contains the harness package classes and variables.
 *
 */
//===========================================================================//

namespace nemesis
{
inline std::string release() { return "x.x.x"; }
}

#endif // harness_Release_hh

//---------------------------------------------------------------------------//
//                        end of harness/Release.hh
//---------------------------------------------------------------------------//
