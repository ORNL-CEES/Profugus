//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Tally.cc
 * \author Thaoms M. Evans
 * \date   Thu May 15 13:27:04 2014
 * \brief  Tally member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Tally.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// TALLY BASE CLASS
//---------------------------------------------------------------------------//
/*!
 * \brief Pure virtual destructor definition.
 */
Tally::~Tally()
{
}

//---------------------------------------------------------------------------//
// TALLY DERIVED CLASSES
//---------------------------------------------------------------------------//
/*!
 * \brief Pure virtual destructor definition.
 */
Source_Tally::~Source_Tally()
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pure virtual destructor definition.
 */
Pathlength_Tally::~Pathlength_Tally()
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pure virtual destructor definition.
 */
Compound_Tally::~Compound_Tally()
{
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Tally.cc
//---------------------------------------------------------------------------//
