//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/TwoGroupCrossSections.hh
 * \author Steven Hamilton
 * \date   Mon Aug 06 14:37:36 2018
 * \brief  TwoGroupCrossSections class declaration.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_TwoGroupCrossSections_hh
#define MC_mc_TwoGroupCrossSections_hh

#include "Assembly_Model.hh"

namespace mc
{

//===========================================================================//
/*!
 * \class TwoGroupCrossSections
 * \brief <++>
 *
 * Long description or discussion goes here.
 */
/*!
 * \example mc/test/tstTwoGroupCrossSections.cc
 *
 * Test of TwoGroupCrossSections.
 */
//===========================================================================//

class TwoGroupCrossSections
{
  public:

    struct XS_Data
    {
        double diffusion[2];
        double absorption[2];
        double scatter;
        double nu_fission[2];
    };

  private:
    // >>> DATA

  public:

    // Constructor
    TwoGroupCrossSections(){}

    // Get XS data at given temperature and density
    XS_Data get_data(Assembly_Model::PIN_TYPE type, double T, double rho);
};

//---------------------------------------------------------------------------//
} // end namespace mc

//---------------------------------------------------------------------------//
#endif // MC_mc_TwoGroupCrossSections_hh

//---------------------------------------------------------------------------//
// end of MC/mc/TwoGroupCrossSections.hh
//---------------------------------------------------------------------------//
