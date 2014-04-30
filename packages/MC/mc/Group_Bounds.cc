//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Group_Bounds.cc
 * \author Thomas M. Evans
 * \date   Wed Apr 30 14:05:49 2014
 * \brief  Group_Bounds member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

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
 * \param p particle type (see mc::Particle_Type).
 */
Group_Bounds::SP_Group_Bounds Group_Bounds::build_logarithmic(double lower,
                                                              double higher,
                                                              int    num_bins)
{
    return std::make_shared<Group_Bounds>(Vec_Dbl());
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
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the upper and lower energy bound of a group.
 */
void Group_Bounds::get_energy(int     group_index,
                              double &lower,
                              double &upper) const
{
    Require(group_index < num_groups());
    Require(group_index + 1 < d_bounds.size());

    upper = d_bounds[group_index];
    lower = d_bounds[group_index + 1];

    Ensure(lower < upper);
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
    Require(energy >= 0.);

    if ((energy > d_bounds.front()) || (energy < d_bounds.back()))
        return false;

    // Find the group index; use std::greater because it's in descending order
    group_index = std::lower_bound(
        d_bounds.begin(), d_bounds.end(), energy, std::greater<double>())
                  - d_bounds.begin() - 1;
    if (group_index == -1)
        ++group_index;
    Check(group_index >= 0);
    Check(group_index < d_bounds.size() - 1);
    return true;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Group_Bounds.cc
//---------------------------------------------------------------------------//
