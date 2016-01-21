//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Constants.hh
 * \author Thomas M. Evans
 * \date   Wed Jul 11 14:17:38 2007
 * \brief  Global constants.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Constants_hh
#define cuda_utils_Constants_hh

#include "math_constants.h"

namespace cuda
{

namespace constants
{

//---------------------------------------------------------------------------//
// BASIC CONSTANTS
//---------------------------------------------------------------------------//

//@{
//! \f$\pi\f$ and \f$\pi\f$-related coefficients.
constexpr double pi          = CUDART_PI;
constexpr double inv_pi      = 1.0 / pi;
constexpr double two_pi      = 2.0 * pi;
constexpr double four_pi     = 4.0 * pi;
constexpr double inv_two_pi  = 1.0 / two_pi;
constexpr double inv_four_pi = 1.0 / four_pi;
//@}

//@{
//! Fractions.
constexpr double one_half    = 1.0 / 2.0;
constexpr double one_third   = 1.0 / 3.0;
constexpr double one_fourth  = 1.0 / 4.0;
constexpr double one_fifth   = 1.0 / 5.0;
constexpr double one_sixth   = 1.0 / 6.0;
constexpr double one_seventh = 1.0 / 7.0;
constexpr double one_eighth  = 1.0 / 8.0;
constexpr double one_tenth   = 1.0 / 10.0;
constexpr double two_thirds  = 2.0 / 3.0;
//@}

//@{
//! Squares and powers.
constexpr double sqrt_two     = CUDART_SQRT_TWO;
//constexpr double sqrt_three   = std::sqrt(3.0);
constexpr double sqrt_four_pi = CUDART_SQRT_2PI * sqrt_two;
constexpr double inv_sqrt_four_pi = 1.0 / sqrt_four_pi;
//@}

//---------------------------------------------------------------------------//
// SCIENCE CONSTANTS
//---------------------------------------------------------------------------//

//! Avagadro's Number (NIST).
constexpr double N_A = 6.0221415e23;
//! Boltzmann Constant (MeV / K) (NIST)
constexpr double k_b = 8.617332478e-11;
//! Mass of a neutron (amu) (NIST)
constexpr double n_mass = 1.0086649160043;
//! Conversion from barns to cm^2
constexpr double barns_per_cm2 = 1.0e24;
//! Conversion from MeV to MJ (NIST).
constexpr double mev2mj = 1.6021766E-19;

} // end namespace cuda::constants

} // end namespace cuda

#endif // cuda_utils_Constants_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Constants.hh
//---------------------------------------------------------------------------//
