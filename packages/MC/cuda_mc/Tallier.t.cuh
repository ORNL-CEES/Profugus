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
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Tallier<Geometry>::Tallier()
{
    d_cell_tally = nullptr;
}

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
void Tallier<Geometry>::set(SDP_Geometry geometry,
                            SDP_Physics  physics)
{
    REQUIRE(geometry.get_host_ptr());
    REQUIRE(geometry.get_device_ptr());
    REQUIRE(physics.get_host_ptr());
    REQUIRE(physics.get_device_ptr());

    d_geometry = geometry.get_device_ptr();
    d_physics  = physics.get_device_ptr();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add a pathlength tally.
 */
template <class Geometry>
void Tallier<Geometry>::add_cell_tally(SDP_Cell_Tally tally)
{
    REQUIRE(tally.get_host_ptr());
    REQUIRE(tally.get_device_ptr());

    // add the tally
    d_cell_tally = tally.get_device_ptr();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tell the tallies to begin active kcode cycles.
 */
template <class Geometry>
void Tallier<Geometry>::begin_active_cycles()
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tell the tallies to begin a new cycle in a kcode calculation.
 */
template <class Geometry>
void Tallier<Geometry>::begin_cycle()
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tell the tallies to end a cycle in a kcode calculation.
 */
template <class Geometry>
void Tallier<Geometry>::end_cycle(double num_particles)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Finalize tallies.
 *
 * Does post-solve processing of tallies including parallel reductions and
 * normalization.
 *
 * \post is_finalized() == true
 */
template <class Geometry>
void Tallier<Geometry>::finalize(double num_particles)
{
    if( d_cell_tally )
        d_cell_tally->finalize(num_particles);
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
void Tallier<Geometry>::reset()
{
    if( d_cell_tally )
        d_cell_tally->reset();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Swap two talliers.
 *
 * This is useful for temporarily deactivating tallying (say, during inactive
 * cycles in a kcode calculation).
 */
template <class Geometry>
void Tallier<Geometry>::swap(Tallier<Geometry> &rhs)
{
    // swap vector internals
    std::swap(d_cell_tally, rhs.d_cell_tally);

    // swap geometry and physics
    std::swap(d_geometry, rhs.d_geometry);
    std::swap(d_physics, rhs.d_physics);
}

} // end namespace cuda_mc

#endif // cuda_mc_Tallier_t_cuh

//---------------------------------------------------------------------------//
//                 end of Tallier.t.cuh
//---------------------------------------------------------------------------//
