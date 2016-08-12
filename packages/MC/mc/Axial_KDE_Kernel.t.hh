//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Axial_KDE_Kernel.t.hh
 * \author Gregory Davidson
 * \date   Mon Feb 16 14:21:15 2015
 * \brief  Axial_KDE_Kernel class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Axial_KDE_Kernel_t_hh
#define MC_mc_Axial_KDE_Kernel_t_hh

#include "Axial_KDE_Kernel.hh"
#include "Sampler.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Geometry>
Axial_KDE_Kernel<Geometry>::Axial_KDE_Kernel(SP_Geometry   geometry,
                                             SP_Physics    physics,
                                             Reject_Method method,
                                             double        coefficient,
                                             double        exponent)
    : Base(geometry, physics, coefficient, exponent)
    , d_method(method)
{   }

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a new position.
 *
 * If the new position is outside the fissionable region, it is rejected.
 */
template<class Geometry>
typename Axial_KDE_Kernel<Geometry>::Space_Vector
Axial_KDE_Kernel<Geometry>::sample_position(const Space_Vector &orig_position,
                                            RNG                &rng) const
{
    REQUIRE(b_physics);
    REQUIRE(b_geometry);
    REQUIRE(rng.assigned());

    if (d_method == FISSION_REJECTION)
    {
        return this->sample_position_fiss_rej(orig_position, rng);
    }
    else if (d_method == CELL_REJECTION)
    {
        return this->sample_position_cell_rej(orig_position, rng);
    }
    else
    {
        throw profugus::assertion("Unknown KDE rejection type");
        // Make compiler happy
        return Space_Vector(0,0,0);
    }
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Sample using fission rejection.
 */
template<class Geometry>
typename Axial_KDE_Kernel<Geometry>::Space_Vector
Axial_KDE_Kernel<Geometry>::sample_position_fiss_rej(
    const Space_Vector &orig_position,
    RNG                &rng) const
{
    REQUIRE(b_physics);
    REQUIRE(b_geometry);
    REQUIRE(rng.assigned());

    // Get the cell at the position
    cell_type cellid = b_geometry->cell(orig_position);

    size_type failures = 0;
    do
    {
        // Sample Epanechnikov kernel
        double epsilon = sampler::sample_epan(rng);

        // Get the bandwidth
        CHECK(b_bndwidth_map.count(cellid) == 1);
        double bandwidth = b_bndwidth_map.find(cellid)->second;
        CHECK(bandwidth >= 0.0);

        // Create a new position
        Space_Vector new_pos(orig_position[def::X],
                             orig_position[def::Y],
                             orig_position[def::Z] + epsilon*bandwidth/2.0);

        // Ensure that the sampled point is in the geometry
        if (b_geometry->boundary_state(new_pos) == geometry::INSIDE)
        {
            // Get the matid for sampled point
            unsigned int matid = b_geometry->matid(new_pos);

            // Get matid from sampled point (may raise error if outside
            // geometry)
            if (b_physics->is_fissionable(b_geometry->matid(new_pos)))
            {
                // Accept: sampled point is fissionable
                b_num_sampled += failures + 1;
                ++b_num_accepted;
                return new_pos;
            }
        }

        // Increment failure counter.
        ++failures;
    } while (failures != 100000);

    // No luck
    b_num_sampled += failures;
    throw profugus::assertion("100000 consecutive nonfissionable rejections in "
                              "KDE [cell " + std::to_string(cellid) + ", "
                              + "original position ("
                              + std::to_string(orig_position[def::X]) + ", "
                              + std::to_string(orig_position[def::Y]) + ", "
                              + std::to_string(orig_position[def::Z]) + ")]");
    return Space_Vector(0,0,0); // Make compiler happy
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample using cell rejection.
 */
template<class Geometry>
typename Axial_KDE_Kernel<Geometry>::Space_Vector
Axial_KDE_Kernel<Geometry>::sample_position_cell_rej(
    const Space_Vector &orig_position,
    RNG                &rng) const
{
    REQUIRE(b_physics);
    REQUIRE(b_geometry);
    REQUIRE(rng.assigned());

    // Get the cell at the position
    cell_type cellid = b_geometry->cell(orig_position);

    size_type failures = 0;
    do
    {
        // Sample Epanechnikov kernel
        double epsilon = sampler::sample_epan(rng);

        // Get the bandwidth if available
        CHECK(b_bndwidth_map.count(cellid) == 1);
        double bandwidth = b_bndwidth_map.find(cellid)->second;
        CHECK(bandwidth >= 0.0);

        // Create a new position
        Space_Vector new_pos(orig_position[def::X],
                             orig_position[def::Y],
                             orig_position[def::Z] + epsilon*bandwidth/2.0);

        // Ensure that the sampled point is in the geometry
        if (b_geometry->boundary_state(new_pos) == geometry::INSIDE)
        {
            // If we are still in the same cell, accept the point, otherwise
            // reject it
            if (b_geometry->cell(new_pos) == cellid)
            {
                // Accept: sampled point is in the same cell
                b_num_sampled += failures + 1;
                ++b_num_accepted;
                return new_pos;
            }
        }

        // Increment failure counter.
        ++failures;
    } while (failures != 100000);

    // No luck
    b_num_sampled += failures;
    throw profugus::assertion("100000 consecutive cell-based rejections in "
                              "KDE kernel [cell " + std::to_string(cellid) +
                              " fission site ("
                              + std::to_string(orig_position[def::X]) + ", "
                              + std::to_string(orig_position[def::Y]) + ", "
                              + std::to_string(orig_position[def::Z]) + ")]");
    return Space_Vector(0,0,0); // Make compiler happy
}

//---------------------------------------------------------------------------//
} // end namespace profugus

#endif // MC_mc_Axial_KDE_Kernel_t_hh

//---------------------------------------------------------------------------//
// end of MC/mc/kde/Axial_KDE_Kernel.t.hh
//---------------------------------------------------------------------------//
