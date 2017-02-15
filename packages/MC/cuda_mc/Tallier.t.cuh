//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Tallier.t.cuh
 * \author Thomas M. Evans
 * \date   Mon May 12 12:15:30 2014
 * \brief  Tallier member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Tallier_t_cuh
#define cuda_mc_Tallier_t_cuh

#include "cuda_utils/CudaDBC.hh"
#include "Tallier.cuh"

namespace cuda_mc
{
//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set the geometry and physics.
 *
 * \param geometry
 * \param physics
 */
template <class Geometry>
void Tallier_DMM<Geometry>::set(SDP_Geometry geometry,
                                SDP_Physics  physics)
{
    d_geometry = geometry;
    d_physics  = physics;
    REQUIRE(d_geometry.get_device_ptr());
    REQUIRE(d_physics.get_device_ptr());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add a cell tally.
 */
template <class Geometry>
void Tallier_DMM<Geometry>::add_cell_tally(SP_Cell_Tally_DMM tally)
{
    REQUIRE(tally);

    // add the tally
    d_cell_tally_dmm = tally;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add a keff tally.
 */
template <class Geometry>
void Tallier_DMM<Geometry>::add_keff_tally(SP_Keff_Tally_DMM tally)
{
    REQUIRE(tally);

    // add the tally
    d_keff_tally_dmm = tally;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tell the tallies to begin active kcode cycles.
 */
template <class Geometry>
void Tallier_DMM<Geometry>::begin_active_cycles()
{
    if (d_keff_tally_dmm)
        d_keff_tally_dmm->begin_active_cycles();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tell the tallies to begin a new cycle in a kcode calculation.
 */
template <class Geometry>
void Tallier_DMM<Geometry>::begin_cycle(def::size_type num_particles)
{
    if (d_keff_tally_dmm)
        d_keff_tally_dmm->begin_cycle(num_particles);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tell the tallies to end a cycle in a kcode calculation.
 */
template <class Geometry>
void Tallier_DMM<Geometry>::end_cycle(double num_particles)
{
    if (d_keff_tally_dmm)
        d_keff_tally_dmm->end_cycle(num_particles);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tell the tallies to end a particle batch
 */
template <class Geometry>
void Tallier_DMM<Geometry>::end_batch(double num_particles)
{
    if (d_cell_tally_dmm)
        d_cell_tally_dmm->end_batch(num_particles);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Finalize tallies.
 *
 * Does post-solve processing of tallies including parallel reductions and
 * normalization.
 */
template <class Geometry>
void Tallier_DMM<Geometry>::finalize(double num_particles)
{
    if (d_cell_tally_dmm)
        d_cell_tally_dmm->finalize(num_particles);

    if (d_keff_tally_dmm)
        d_keff_tally_dmm->finalize(num_particles);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reset results in all tallies
 *
 * This does not remove tallies from the tallier: it will just reset the
 * accumulators. See implementation of the Tally daughter classes for details,
 * but generally this doesn't clear the values in existing smart-pointer
 * fields.
 */
template <class Geometry>
void Tallier_DMM<Geometry>::reset()
{
    if (d_cell_tally_dmm)
        d_cell_tally_dmm->reset();

    if (d_keff_tally_dmm)
        d_keff_tally_dmm->reset();
}

} // end namespace cuda_mc

#endif // cuda_mc_Tallier_t_cuh

//---------------------------------------------------------------------------//
//                 end of Tallier.t.cuh
//---------------------------------------------------------------------------//
