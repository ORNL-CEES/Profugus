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
    COLLISION = 0,          //!< Collided
    BOUNDARY = 1,           //!< Hit internal boundary
    SCATTER = 2,            //!< Scattered
    IMPLICIT_CAPTURE = 3,   //!< Weight decreased through implicit absorption
    WEIGHT_WINDOW = 4,      //!< Passed through weight window successfully
    ROULETTE_SURVIVE = 5,   //!< Survived Russian roulette
    BOUNDARY_MESH = 6,      //!< Encountered domain decomposition boundary
    BORN = 7,               //!< Born
    SPLIT = 8,              //!< Split by weight window
                       
    STILL_ALIVE = 9,        //!< Particle is still alive. Used for event-based sort.

    CUTOFF = 10,            //!< Cutoff by energy
    ABSORPTION = 11,        //!< Was absorbed
    ESCAPE = 12,            //!< Left problem through external boundary
    ERROR_KILLED = 13,      //!< Encountered unexpected error
    ROULETTE_KILLED = 14,   //!< Killed by Russian roulette
    END_EVENT = 15
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
