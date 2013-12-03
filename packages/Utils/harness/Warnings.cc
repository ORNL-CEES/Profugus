//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/Warnings.cc
 * \author Thomas M. Evans
 * \date   Sun Feb 26 20:54:46 2012
 * \brief  Warnings member definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Warnings.hh"
#include "DBC.hh"

namespace nemesis
{

//---------------------------------------------------------------------------//
// EXTERNAL DEFINITIONS OF WARNINGS CLASS
//---------------------------------------------------------------------------//

namespace warn
{

Warnings warnings;

} // end of warn

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
Warnings::Warnings()
{
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Add a warning.
 */
void Warnings::add(const std::string &w)
{
    d_warnings.push_back(w);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Access the oldest warning.
 *
 * This can only be called if \c empty() is \c false .
 */
std::string Warnings::pop()
{
    Require(!d_warnings.empty());

    std::string w(d_warnings.front());
    d_warnings.pop_front();

    return w;
}

} // end namespace nemesis

//---------------------------------------------------------------------------//
//                 end of Warnings.cc
//---------------------------------------------------------------------------//
