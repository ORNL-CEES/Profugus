//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Definitions.hh
 * \author Thomas M. Evans
 * \date   Friday April 25 16:46:37 2014
 * \brief  Monte Carlo Definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Definitions_hh
#define mc_Definitions_hh

namespace profugus
{

//---------------------------------------------------------------------------//

namespace events
{

//! Monte Carlo event descriptors.
enum Event {
    COLLISION = 0,    //!< Collided
    BOUNDARY,         //!< Hit internal boundary
    SCATTER,          //!< Scattered
    IMPLICIT_CAPTURE, //!< Weight decreased through implicit absorption
    WEIGHT_WINDOW,    //!< Passed through weight window successfully
    ROULETTE_SURVIVE, //!< Survived Russian roulette
    BOUNDARY_MESH,    //!< Encountered domain decomposition boundary
    BORN,             //!< Born
    SPLIT,            //!< Split by weight window
                       
    STILL_ALIVE,      //!< Particle is still alive. Used for event-based sort.

    CUTOFF,           //!< Cutoff by energy
    ABSORPTION,       //!< Was absorbed
    ESCAPE,           //!< Left problem through external boundary
    ERROR_KILLED,     //!< Encountered unexpected error
    ROULETTE_KILLED,  //!< Killed by Russian roulette
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

#endif // mc_Definitions_hh

//---------------------------------------------------------------------------//
//              end of mc/Definitions.hh
//---------------------------------------------------------------------------//
