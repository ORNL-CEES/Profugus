//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/kde/KDE_Kernel_Resample.cc
 * \author Gregory Davidson
 * \date   Mon Feb 16 14:21:15 2015
 * \brief  KDE_Kernel_Resample class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "KDE_Kernel_Resample.hh"
#include "Sampler.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
KDE_Kernel_Resample::KDE_Kernel_Resample(SP_Geometry geometry,
                                         SP_Physics  physics,
                                         double      coefficient,
                                         double      exponent)
    : Base(geometry, physics, coefficient, exponent)
{   }

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a new position.
 *
 * If the new position is outside the fissionable region, it is rejected.
 */
KDE_Kernel_Resample::Space_Vector
KDE_Kernel_Resample::sample_position(const Space_Vector &orig_position,
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
    } while (failures != 1000);

    // No luck
    b_num_sampled += failures;
    std::cout << "1000 consecutive nonfissionable rejections in KDE."
              << std::endl;
    return Space_Vector(0,0,0);
}

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
// end of MC/mc/kde/KDE_Kernel_Resample.cc
//---------------------------------------------------------------------------//
