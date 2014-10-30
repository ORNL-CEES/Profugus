//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   acc/Physics.acc.cc
 * \author Seth R Johnson
 * \date   Wed Oct 29 10:32:32 2014
 * \brief  Physics class definitions.
 * \note   Copyright (c) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Physics.hh"

#include "harness/DBC.hh"

#include "Change_Direction.hh"

namespace acc
{
//---------------------------------------------------------------------------//
void Physics::complete()
{
    REQUIRE(d_num_mats > 0);
    REQUIRE(d_num_groups > 0);
    REQUIRE(d_total != 0);
    REQUIRE(d_nusigf != 0);
    REQUIRE(d_scatter != 0);
    REQUIRE(d_scatter_ratio != 0);
    REQUIRE(d_fissionable != 0);
    int ne = dv_total.size();
    int ns = dv_scatter.size();
#pragma acc enter data \
    copyin(this)
#pragma acc enter data \
    copyin(d_total[0:ne], d_nusigf[0:ne], d_scatter[0:ns],\
           d_scatter_ratio[0:ne], d_fissionable[0:d_num_mats])
}

//---------------------------------------------------------------------------//
Physics::~Physics()
{
#pragma acc exit data \
    delete(d_total, d_nusigf, d_scatter, d_scatter_ratio, d_fissionable)
#pragma acc exit data \
    delete(this)
}

//---------------------------------------------------------------------------//
// PHYSICS
//---------------------------------------------------------------------------//
/*!
 * \brief Total macro XS
 */
double Physics::total(int matid, int group) const
{
    return d_total[vector_index(matid, group)];
}

//---------------------------------------------------------------------------//
//! Scattering ratio
double Physics::scattering_ratio(int matid, int group) const
{
    return d_scatter_ratio[vector_index(matid, group)];
}

//---------------------------------------------------------------------------//
//! Nu fission
double Physics::nusigf(int matid, int group) const
{
    return d_nusigf[vector_index(matid, group)];
}

//---------------------------------------------------------------------------//
void Physics::collide(Particle& p) const
{
    const int matid = p.matid;
    const int group = p.group;

    // Get scattering ratio
    const double c = scattering_ratio(matid, group);

    // Adjust particle weight by scattering ratio; if it's a pure absorber,
    // this will set to zero
    p.wt *= c;

    // >>> Sample exit group
    const double* outscatter = d_scatter + matrix_index(matid, 0, group);
    //CHECK(outscatter + 1 == d_scatter + matrix_index(matid, 1, group));

    const double xi = 0.5;
    const double threshold = xi * total(matid, group) * c;
    double accum = 0;

    // Build CDF inline
    int gp = 0;
    for (int gp_end = d_num_groups - 1; gp != gp_end; ++gp)
    {
        accum += *outscatter;
        if (accum >= threshold)
            break;

        ++outscatter;
    }

    // Set the new group
    p.group = gp;

    // >>> Scatter
    {
        double costheta = 0.5;
        double phi = 0.01;
        change_direction(costheta, phi, p.geo_state);
    }
}

//---------------------------------------------------------------------------//
//! Index into vector data
int Physics::vector_index(int mat, int group) const
{
    return mat * d_num_groups + group;
}

//---------------------------------------------------------------------------//
//! Index into matrix data
int Physics::matrix_index(int mat, int out_group, int in_group) const
{
    return (mat * d_num_groups + in_group) * d_num_groups + out_group;
}

//---------------------------------------------------------------------------//
} // end namespace acc

//---------------------------------------------------------------------------//
// end of acc/Physics.acc.cc
//---------------------------------------------------------------------------//
