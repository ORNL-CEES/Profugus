//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   driver/Manager.hh
 * \author Thomas M. Evans
 * \date   Fri Mar 14 11:32:36 2014
 * \brief  Manager class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef driver_Manager_hh
#define driver_Manager_hh

#include "Problem_Builder.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Manager
 * \brief Manager class that drives the SPN miniapp.
 */
//===========================================================================//

class Manager
{
  private:
    // >>> DATA

    // Problem builder.
    Problem_Builder d_builder;

  public:
    // Constructor.
    Manager();

    // Setup the problem.
    void setup(const std::string &xml_file);
};

} // end namespace profugus

#endif // driver_Manager_hh

//---------------------------------------------------------------------------//
//                 end of Manager.hh
//---------------------------------------------------------------------------//
