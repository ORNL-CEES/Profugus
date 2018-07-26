//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc_driver/Water_Properties.hh
 * \author Steven Hamilton
 * \date   Wed Jul 25 14:29:11 2018
 * \brief  Water_Properties class declaration.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_driver_Water_Properties_hh
#define MC_mc_driver_Water_Properties_hh

namespace mc
{

//===========================================================================//
/*!
 * \class Water_Properties
 * \brief Get material properties for water
 */
/*!
 * \example mc_driver/test/tstWater_Properties.cc
 *
 * Test of Water_Properties.
 */
//===========================================================================//

class Water_Properties
{
  public:

    static double Density(double h, double p);
    static double Enthalpy(double T, double p);
    static double Temperature(double h, double p);
};

//---------------------------------------------------------------------------//
} // end namespace mc

//---------------------------------------------------------------------------//
#endif // MC_mc_driver_Water_Properties_hh

//---------------------------------------------------------------------------//
// end of MC/mc_driver/Water_Properties.hh
//---------------------------------------------------------------------------//
