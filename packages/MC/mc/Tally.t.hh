//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Tally.t.hh
 * \author Thaoms M. Evans
 * \date   Thu May 15 13:27:04 2014
 * \brief  Tally template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Tally_t_hh
#define mc_Tally_t_hh

#include "Tally.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// TALLY BASE CLASS
//---------------------------------------------------------------------------//
/*!
 * \brief Pure virtual destructor definition.
 */
template <class Geometry>
Tally<Geometry>::~Tally()
{
}

//---------------------------------------------------------------------------//
// TALLY DERIVED CLASSES
//---------------------------------------------------------------------------//
/*!
 * \brief Pure virtual destructor definition.
 */
template <class Geometry>
Source_Tally<Geometry>::~Source_Tally()
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pure virtual destructor definition.
 */
template <class Geometry>
Pathlength_Tally<Geometry>::~Pathlength_Tally()
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pure virtual destructor definition.
 */
template <class Geometry>
Compound_Tally<Geometry>::~Compound_Tally()
{
}

} // end namespace profugus

#endif // mc_Tally_t_hh

//---------------------------------------------------------------------------//
//                 end of Tally.t.hh
//---------------------------------------------------------------------------//
