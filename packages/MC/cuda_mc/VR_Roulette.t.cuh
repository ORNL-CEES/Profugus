//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/VR_Roulette.t.cuh
 * \author Thomas M. Evans
 * \date   Fri May 09 13:09:37 2014
 * \brief  VR_Roulette member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_VR_Roulette_t_cuh
#define cuda_mc_VR_Roulette_t_cuh

#include "VR_Roulette.cuh"
#include "mc/Definitions.hh"
#include "harness/DBC.hh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
VR_Roulette<Geometry>::VR_Roulette(RCP_Std_DB db)
{
    REQUIRE(!db.is_null());

    // Get default weight cutoffs and survival
    d_Wc = db->get("weight_cutoff", 0.25);
    d_Ws = db->get("weight_survival", 2.0 * d_Wc);

    INSIST(d_Wc >= 0, "Weight cutoff must be nonnegative");

    INSIST(d_Ws > 0 ? d_Ws > d_Wc : d_Wc == 0.,
            "Survival weight must be greater than cutoff weight");
}

} // end namespace cuda_mc 

#endif // cuda_mc_VR_Roulette_t_cuh

//---------------------------------------------------------------------------//
//                 end of VR_Roulette.t.cuh
//---------------------------------------------------------------------------//
