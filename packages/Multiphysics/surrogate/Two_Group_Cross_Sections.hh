//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Two_Group_Cross_Sections.hh
 * \author Steven Hamilton
 * \date   Mon Aug 06 14:37:36 2018
 * \brief  Two_Group_Cross_Sections class declaration.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Two_Group_Cross_Sections_hh
#define MC_mc_Two_Group_Cross_Sections_hh

#include "Assembly_Model.hh"

namespace mc
{

//===========================================================================//
/*!
 * \class Two_Group_Cross_Sections
 * \brief Class for interpolated two group cross section data.
 */
/*!
 * \example mc/test/tstTwo_Group_Cross_Sections.cc
 *
 * Test of Two_Group_Cross_Sections.
 */
//===========================================================================//

class Two_Group_Cross_Sections
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
    Two_Group_Cross_Sections(){}

    // Get XS data at given temperature and density
    XS_Data get_data(Assembly_Model::PIN_TYPE type, double T, double rho);
};

//---------------------------------------------------------------------------//
} // end namespace mc

//---------------------------------------------------------------------------//
#endif // MC_mc_Two_Group_Cross_Sections_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Two_Group_Cross_Sections.hh
//---------------------------------------------------------------------------//
