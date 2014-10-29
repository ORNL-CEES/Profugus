//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   acc/Physics.cpu.cc
 * \author Seth R Johnson
 * \date   Wed Oct 29 11:07:47 2014
 * \brief  Physics.cpu class definitions.
 * \note   Copyright (c) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Physics.hh"

#include "xs/XS.hh"

using namespace profugus;

namespace acc
{
//---------------------------------------------------------------------------//
Physics::Physics(const XS& xs)
{
    d_num_mats = xs.num_mat();
    d_num_groups = xs.num_groups();
    dv_total.resize(num_vector_elements());
    dv_nusigf.resize(num_vector_elements());
    dv_scatter.resize(num_matrix_elements());

    // Copy data for every material into our flattened structures
    for (int m = 0; m < d_num_mats; ++m)
    {
        const XS::Matrix& scatter = xs.matrix(m, 0); // Pn order zero
        const XS::Vector& total   = xs.vector(m, XS::TOTAL);
        const XS::Vector& nusigf  = xs.vector(m, XS::NU_SIG_F);

        for (int g = 0; g < d_num_groups; ++g)
        {
            const int veci = vector_index(m, g);
            CHECK(veci < dv_total.size());
            CHECK(veci < dv_nusigf.size());

            // Copy vector data
            dv_total[veci]  = total[g];
            dv_nusigf[veci] = nusigf[g];

            // Loop over outscattering groups
            for (int gp = 0; gp < d_num_groups; ++gp)
            {
                // Index in scattering from g -> g'
                const int mati = matrix_index(m, gp, g);
                CHECK(mati < dv_scatter.size());

                dv_scatter[mati] = scatter(gp, g);
            }
        }
    }

    // Set pointers
    d_total = dv_total.data();
    d_nusigf = dv_nusigf.data();
    d_scatter = dv_scatter.data();
}

//---------------------------------------------------------------------------//
} // end namespace acc

//---------------------------------------------------------------------------//
// end of acc/Physics.cpu.cc
//---------------------------------------------------------------------------//
