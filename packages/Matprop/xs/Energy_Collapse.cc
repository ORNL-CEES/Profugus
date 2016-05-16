//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Matprop/xs/Energy_Collapse.cc
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 13:39:38 2014
 * \brief  Energy_Collapse member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <numeric>
#include "Energy_Collapse.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Collapse material database from fine to coarse.
 */
Energy_Collapse::RCP_Mat_DB Energy_Collapse::collapse_all_mats(
    RCP_Mat_DB     fine_mat,
    const Vec_Int &collapse_vec,
    const Vec_Dbl &weights)
{
    typedef Vec_Int::const_iterator mat_iter;
    typedef Mat_DB_t::XS_t          XS;

    REQUIRE( fine_mat->xs().num_groups() == weights.size() );
    REQUIRE( fine_mat->xs().num_groups() ==
             std::accumulate(collapse_vec.begin(),collapse_vec.end(),0) );
    REQUIRE( fine_mat->xs().num_groups() >= 2 );

    int fine_grps   = fine_mat->xs().num_groups();
    int coarse_grps = collapse_vec.size();

    // Get starting fine group indices for each coarse group
    Vec_Int start_ind(coarse_grps);
    start_ind[0] = 0;
    if (coarse_grps > 1)
    {
        std::partial_sum(collapse_vec.begin(), collapse_vec.end()-1,
                         start_ind.begin()+1);
    }

    // get the cross sections from the fine_mat
    const XS &xs = fine_mat->xs();

    // get the matids in the database
    Vec_Int matids;
    xs.get_matids(matids);

    // material id
    int m = 0;

    // pn order
    int pn_order = xs.pn_order();

    // create coarse group cross sections
    Mat_DB_t::RCP_XS xsc = Teuchos::rcp(new XS);
    xsc->set(pn_order, coarse_grps);

    // coarse group totals
    XS::OneDArray totc(coarse_grps, 0.0);
    XS::TwoDArray sctc(coarse_grps, coarse_grps, 0.0);

    // Process each material in original
    for (mat_iter mitr = matids.begin(); mitr != matids.end(); ++mitr)
    {
        // set the material id
        m = *mitr;

        // get the total cross sections for the fine group
        const XS::Vector &totf = xs.vector(m, XS::TOTAL);

        // loop over coarse groups
        for (int gc = 0; gc < coarse_grps; ++gc)
        {
            // Get fine first and last for this coarse
            int g_first = start_ind[gc];
            int g_last  = g_first + collapse_vec[gc];

            // Process total cross section
            double sigma = 0.0;
            double g_sum = 0.0;
            for (int gf = g_first; gf < g_last; ++gf)
            {
                sigma += totf[gf] * weights[gf];
                g_sum += weights[gf];
            }
            sigma    /= g_sum;
            totc[gc]  = sigma;
        }

        // add the coarse total cross sections
        xsc->add(m, XS::TOTAL, totc);

        // add the scattering cross sections
        for (int l = 0; l <= pn_order; ++l)
        {
            // get the scattering cross section for the fine group
            const XS::Matrix &sctf = xs.matrix(m, l);

            for (int gc = 0; gc < coarse_grps; ++gc)
            {
                // Get fine first and last for this coarse
                int g_first = start_ind[gc];
                int g_last  = g_first + collapse_vec[gc];

                for (int gpc = 0; gpc < coarse_grps; ++gpc)
                {
                    int gp_first  = start_ind[gpc];
                    int gp_last   = gp_first + collapse_vec[gpc];
                    double sig_s  = 0.0;
                    double gp_sum = 0.0;
                    for (int gpf = gp_first; gpf < gp_last; ++gpf)
                    {
                        gp_sum += weights[gpf];
                        for (int gf = g_first; gf < g_last; ++gf)
                        {
                            sig_s += sctf(gf, gpf) * weights[gpf];
                        }
                    }
                    sig_s         /= gp_sum;
                    sctc(gc, gpc)  = sig_s;
                }
            }

            // add the coarse scattering cross sections
            xsc->add(m, l, sctc);
        }
    }

    // complete the coarse cross sections
    xsc->complete();

    // make the coarse mat-db
    RCP_Mat_DB coarse_mat = Teuchos::rcp(new Mat_DB_t);
    coarse_mat->set(xsc, fine_mat->num_cells());

    // Assign matids to cells
    for (int cell = 0, Nc = fine_mat->num_cells(); cell < Nc; ++cell)
    {
        coarse_mat->matid(cell) = fine_mat->matid(cell);
    }

    return coarse_mat;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Energy_Collapse.cc
//---------------------------------------------------------------------------//
