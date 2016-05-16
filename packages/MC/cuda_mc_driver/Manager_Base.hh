//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc_driver/Manager_Base.hh
 * \author Steven Hamilton
 * \date   Wed Nov 25 11:20:35 2015
 * \brief  Manager_Base class definition.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_driver_Manager_Base_hh
#define cuda_mc_driver_Manager_Base_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Manager_Base
 * \brief Base class for templated managers
 */
//===========================================================================//

class Manager_Base
{
  public:

      typedef Teuchos::RCP<Teuchos::ParameterList> RCP_ParameterList;

      // Constructor
      Manager_Base(){};

      // Virtual Destructor
      virtual ~Manager_Base(){};

      // Setup the problem.
      virtual void setup(RCP_ParameterList master) = 0;

      // Solve the problem.
      virtual void solve() = 0;

      // Output data.
      virtual void output() = 0;
};

} // end namespace cuda_mc

#endif // cuda_mc_driver_Manager_Base_hh

//---------------------------------------------------------------------------//
//                 end of Manager_Base.hh
//---------------------------------------------------------------------------//
