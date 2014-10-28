//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/mc/Group_Bounds.cc
 * \author Thomas M. Evans
 * \date   Wed Apr 30 14:05:49 2014
 * \brief  Group_Bounds member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <algorithm>

#include "Group_Bounds.hh"

#include "harness/DBC.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// STATIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Construct a logarithmic group bounds structure.
 *
 * \param lower lower energy bound (eV).
 * \param higher higher energy bound (eV).
 * \param num_bins number of logarithmic energy bins.
 */
Group_Bounds::SP_Group_Bounds Group_Bounds::build_logarithmic(double lower,
                                                              double higher,
                                                              int    num_bins)
{
    INSIST(lower > 0.0, "Lower energy bound must be > ZERO");
    INSIST(higher > lower,
            "Upper energy bound must be greater than the lower energy bound");
    INSIST(num_bins > 0, "Must have a least one energy bin");

    // allocate the energy bounds
    Vec_Dbl energy_mesh(num_bins + 1, 0.0);

    // calculate the width of each energy bin
    double de = std::log(higher / lower) / num_bins;

    // assign the energy bounds
    energy_mesh[0]        = higher;
    energy_mesh[num_bins] = lower;
    for (int g = 1; g < num_bins; ++g)
        energy_mesh[g] = higher / std::exp(g * de);

    // make the group bounds
    SP_Group_Bounds gb(std::make_shared<Group_Bounds>(energy_mesh));

    ENSURE(gb);
    ENSURE(gb->num_groups() == num_bins);
    return gb;
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Construct with neutron and photon boundaries.
 *
 * Boundaries must be monotonic decreasing or empty. At least one group
 * boundaries must be provided.
 */
Group_Bounds::Group_Bounds(const Vec_Dbl &bounds)
    : d_bounds(bounds)
{
    // check monotonicity
    for (int g = 0; g < num_groups(); ++g)
    {
        VALIDATE(d_bounds[g] > d_bounds[g+1],
                  "Energy bounds for group " << g << " <= bounds at "
                  << g+1);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the upper and lower energy bound of a group.
 */
void Group_Bounds::get_energy(int     group_index,
                              double &lower,
                              double &upper) const
{
    REQUIRE(group_index < num_groups());
    REQUIRE(group_index + 1 < d_bounds.size());

    upper = d_bounds[group_index];
    lower = d_bounds[group_index + 1];

    ENSURE(lower < upper);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the flattened group index for a particle at a given energy
 *
 * Lower bounds of an energy group are included as part of that group.
 *
 * \param[in] energy       Particle energy
 * \param[out] group_index Group index in multigroup data
 *
 * \return true if found, false if not
 */
bool Group_Bounds::find(const double  energy,
                        int          &group_index) const
{
    REQUIRE(energy >= 0.);

    if ((energy > d_bounds.front()) || (energy < d_bounds.back()))
        return false;

    // Find the group index; use std::greater because it's in descending order
    group_index = std::lower_bound(
        d_bounds.begin(), d_bounds.end(), energy, std::greater<double>())
                  - d_bounds.begin() - 1;
    if (group_index == -1)
        ++group_index;
    CHECK(group_index >= 0);
    CHECK(group_index < d_bounds.size() - 1);
    return true;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Group_Bounds.cc
//---------------------------------------------------------------------------//
