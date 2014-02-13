//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Constants.hh
 * \author Thomas M. Evans
 * \date   Wed Jul 11 14:17:38 2007
 * \brief  Global constants.
 * \note   Copyright (C) 2007 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Constants_hh
#define utils_Constants_hh

#include <cmath>
#include <limits>

namespace profugus
{

namespace constants
{

//---------------------------------------------------------------------------//
// BASIC CONSTANTS
//---------------------------------------------------------------------------//

//@{
//! \f$\pi\f$ and \f$\pi\f$-related coefficients.
const double pi          = 2.0 * std::asin(1.0);
const double inv_pi      = 1.0 / pi;
const double two_pi      = 2.0 * pi;
const double four_pi     = 4.0 * pi;
const double inv_two_pi  = 1.0 / two_pi;
const double inv_four_pi = 1.0 / four_pi;
//@}

//@{
//! Fractions.
const double one_half    = 1.0 / 2.0;
const double one_third   = 1.0 / 3.0;
const double one_fourth  = 1.0 / 4.0;
const double one_fifth   = 1.0 / 5.0;
const double one_sixth   = 1.0 / 6.0;
const double one_seventh = 1.0 / 7.0;
const double one_eighth  = 1.0 / 8.0;
const double one_tenth   = 1.0 / 10.0;
const double two_thirds  = 2.0 / 3.0;
//@}

//@{
//! Squares and powers.
const double sqrt_two     = std::sqrt(2.0);
const double sqrt_three   = std::sqrt(3.0);
const double sqrt_four_pi = std::sqrt(four_pi);
const double inv_sqrt_four_pi = 1.0 / sqrt_four_pi;
//@}

//@{
//! Numeric limits.
const double tiny = std::numeric_limits<double>::min();
const double huge = std::numeric_limits<double>::max();
//@}

//---------------------------------------------------------------------------//
// SCIENCE CONSTANTS
//---------------------------------------------------------------------------//

//! Avagadro's Number (NIST).
const double N_A = 6.0221415e23;
//! Boltzmann Constant (MeV / K) (NIST)
const double k_b = 8.617332478e-11;
//! Mass of a neutron (amu) (NIST)
const double n_mass = 1.0086649160043;
//! Conversion from barns to cm^2
const double barns_per_cm2 = 1.0e24;
//! Conversion from MeV to MJ (NIST).
const double mev2mj = 1.6021766E-19;

} // end namespace profugus::constants

} // end namespace profugus

#endif // utils_Constants_hh

//---------------------------------------------------------------------------//
//              end of utils/Constants.hh
//---------------------------------------------------------------------------//
