//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   driver/Manager.cc
 * \author Thomas M. Evans
 * \date   Fri Mar 14 11:32:36 2014
 * \brief  Manager member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "harness/DBC.hh"
#include "Manager.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Manager::Manager()
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Setup the problem.
 * \param xml_file
 */
void Manager::setup(const std::string &xml_file)
{
    // use the problem builder to setup the problem
    d_builder.setup(xml_file);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Manager.cc
//---------------------------------------------------------------------------//
