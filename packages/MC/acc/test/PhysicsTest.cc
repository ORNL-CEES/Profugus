//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   acc/test/PhysicsTest.cc
 * \author Seth R Johnson
 * \date   Thu Oct 30 09:35:40 2014
 * \brief  PhysicsTest class definitions.
 * \note   Copyright (c) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "PhysicsTest.hh"

#include "../RNG.hh"

void calc_total(
        const acc::Physics& physics,
        const int*     matids,
        const int*     groups,
        int            num_particles,
        double*        total)
{
#pragma acc parallel loop \
        present(physics) \
        copyin(groups[0:num_particles], matids[0:num_particles]) \
        copyout(total[0:num_particles])
    for (int i = 0; i < num_particles; ++i)
    {
        total[i] = physics.total(matids[i], groups[i]);
    }
}

// With start groups and matids, loop over particles and accumulate weights
void multi_scatter(
        const acc::Physics&  physics,
        const int*           matids,
        const int*           start_groups,
        int                  num_particles,
        int                  num_steps,
        int*                 end_group,
        double*              end_wt)
{
//#pragma acc parallel loop \
//        present(physics) \
//        copyin(start_groups[0:num_particles], matids[0:num_particles]) \
//        pcopyout(end_group[0:num_particles], end_wt[0:num_particles])
    for (int n = 0; n < num_particles; ++n)
    {
        // Initialize particle on GPU
        acc::Particle p;
        p.matid = matids[n];
        p.group = start_groups[n];
        p.wt = 1;
        p.geo_state.dir[0] = p.geo_state.dir[1] = 0;
        p.geo_state.dir[2] = 1;

        // Loop over steps
#pragma acc loop seq
        for (int s = 0; s < num_steps; ++s)
        {
            physics.collide(p);
        }
        end_group[n] = p.group;
        end_wt[n]    = p.wt;
    }

}

//---------------------------------------------------------------------------//
// end of acc/test/PhysicsTest.cc
//---------------------------------------------------------------------------//
