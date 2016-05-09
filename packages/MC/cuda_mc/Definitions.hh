//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Definitions.hh
 * \author Thomas M. Evans
 * \date   Friday April 25 16:46:37 2014
 * \brief  Monte Carlo Definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Definitions_hh
#define cuda_mc_Definitions_hh

namespace cuda_profugus
{

//---------------------------------------------------------------------------//

namespace events
{

//! Monte Carlo event descriptors.
enum Event {
    COLLISION = 0,		//!< Collided
    ABSORPTION,			//!< Was absorbed
    SCATTER,			//!< Scattered
    BOUNDARY,			//!< Hit internal boundary
    CUTOFF,			//!< Cutoff by energy
    ESCAPE,			//!< Left problem through external boundary
    ERROR_KILLED,		//!< Encountered unexpected error
    IMPLICIT_CAPTURE,		//!< Weight decreased through implicit absorption
    VR_POST_SURFACE,		//!< Ready for post-surface variance reduction
    ROULETTE_KILLED,		//!< Killed by Russian roulette
    ROULETTE_SURVIVE,		//!< Survived Russian roulette
    SPLIT,			//!< Split by weight window
    WEIGHT_WINDOW,		//!< Passed through weight window successfully
    BOUNDARY_MESH,		//!< Encountered domain decomposition boundary
    TAKE_STEP,                  //!< Take a transport step
    BORN,			//!< Born
    DEAD,			//!< The particle is dead and has no more events
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

} // end namespace profugus

#endif // cuda_mc_Definitions_hh

//---------------------------------------------------------------------------//
//              end of cuda_mc/Definitions.hh
//---------------------------------------------------------------------------//
