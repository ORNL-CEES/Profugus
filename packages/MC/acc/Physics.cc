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
    dv_fissionable.resize(d_num_mats);
    dv_total.resize(num_vector_elements());
    dv_nusigf.resize(num_vector_elements());
    dv_scatter_ratio.resize(num_vector_elements());
    dv_outscatter_pdf.resize(num_matrix_elements());

    // Copy data for every material into our flattened structures
    for (int m = 0; m < d_num_mats; ++m)
    {
        const XS::Matrix& scatter = xs.matrix(m, 0); // Pn order zero
        const XS::Vector& total   = xs.vector(m, XS::TOTAL);
        const XS::Vector& nusigf  = xs.vector(m, XS::NU_SIG_F);

        bool is_fiss = false;
        for (int g = 0; g < d_num_groups; ++g)
        {
            const int veci = vector_index(m, g);
            CHECK(veci < dv_total.size());
            CHECK(veci < dv_nusigf.size());

            // Copy vector data
            dv_total[veci]  = total[g];
            dv_nusigf[veci] = nusigf[g];

            // Check for fissionable
            if (nusigf[g] > 0.0)
               is_fiss = true;
            // Loop over outscattering groups
            double tot_scat = 0;
            for (int gp = 0; gp < d_num_groups; ++gp)
            {
                tot_scat += scatter(gp, g);
            }
            for (int gp = 0; gp < d_num_groups; ++gp)
            {
                // Index in scattering from g -> g'
                const int mati = matrix_index(m, gp, g);
                CHECK(mati < dv_outscatter_pdf.size());

                dv_outscatter_pdf[mati] = scatter(gp, g) / tot_scat;
            }
            // Compute and set scattering ratio
            dv_scatter_ratio[veci] = tot_scat / total[g];
        }
        // Set fissionable flag for this material
        dv_fissionable[m] = is_fiss;
    }

    // Set pointers
    d_total = dv_total.data();
    d_nusigf = dv_nusigf.data();
    d_outscatter_pdf = dv_outscatter_pdf.data();
    d_scatter_ratio = dv_scatter_ratio.data();
    d_fissionable = dv_fissionable.data();

    // Complete to copy data to GPU
    complete();
}

//---------------------------------------------------------------------------//
} // end namespace acc

//---------------------------------------------------------------------------//
// end of acc/Physics.cpu.cc
//---------------------------------------------------------------------------//
