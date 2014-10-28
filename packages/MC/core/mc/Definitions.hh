//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/mc/Definitions.hh
 * \author Thomas M. Evans
 * \date   Friday April 25 16:46:37 2014
 * \brief  Monte Carlo Definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef core_mc_Definitions_hh
#define core_mc_Definitions_hh

namespace profugus
{

//---------------------------------------------------------------------------//

namespace events
{

//! Monte Carlo event descriptors.
enum Event {
    COLLISION = 0,    //!< Collided
    ABSORPTION,       //!< Was absorbed
    SCATTER,          //!< Scattered
    BOUNDARY,         //!< Hit internal boundary
    CUTOFF,           //!< Cutoff by energy
    ESCAPE,           //!< Left problem through external boundary
    ERROR_KILLED,     //!< Encountered unexpected error
    IMPLICIT_CAPTURE, //!< Weight decreased through implicit absorption
    ROULETTE_KILLED,  //!< Killed by Russian roulette
    ROULETTE_SURVIVE, //!< Survived Russian roulette
    SPLIT,            //!< Split by weight window
    WEIGHT_WINDOW,    //!< Passed through weight window successfully
    BOUNDARY_MESH,    //!< Encountered domain decomposition boundary
    BORN,             //!< Born
    END_EVENT
};

} // end namespace events

//---------------------------------------------------------------------------//

namespace physics
{

//! Physics reaction types.
enum Reaction_Type {
    FLUX       = -1, //!< i.e. no reaction type
    TOTAL      = 0,  //!< Total of all reactions
    ABSORPTION,      //!< Absorption reaction
    SCATTERING,      //!< Elastic scattering reaction
    FISSION,         //!< Fission reaction
    NU_FISSION,      //!< Average neutron production
    KAPPA_SIGMA,     //!< Recoverable energy production
    N_DIS,           //!< Neutron disappearance
    RAD_CAPTURE,     //!< Radiative capture reaction
    END_REACTION_TYPE
};

} // end namespace physics

} // end namespace mc

#endif // core_mc_Definitions_hh

//---------------------------------------------------------------------------//
//              end of mc/Definitions.hh
//---------------------------------------------------------------------------//
