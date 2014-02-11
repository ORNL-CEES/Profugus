//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Moment_Coefficients.cc
 * \author Thomas M. Evans
 * \date   Mon Feb 10 13:07:08 2014
 * \brief  Moment_Coefficients member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <set>
#include <cmath>

#include "harness/DBC.hh"
#include "comm/global.hh"
#include "comm/Timing.hh"
#include "utils/Definitions.hh"
#include "Moment_Coefficients.hh"

// Anonymous namespace function to convert a moment/matid pair to a single
// size_type for hashing
namespace
{
    def::size_type to_size_type(int mom, int mat)
    {
        return mom + 8 * mat;
    }
}

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Moment_Coefficients::Moment_Coefficients(RCP_ParameterList db,
                                         RCP_Dimensions    dim,
                                         RCP_Mat_DB        mat)


    : d_dim(dim)
    , d_mat(mat)
    , d_Ng(mat->xs().num_groups())
    , d_W(d_Ng, d_Ng)
    , d_c(4, std::vector<Coefficients>(4))
    , d_outscatter_correction(false)
    , d_work(d_Ng)
    , d_ipiv(d_Ng)
{
    Require (!d_dim.is_null());
    Require (!d_mat.is_null());
    Require (d_c.size() == 4);
    Require (d_c[0].size() == 4);
    Require (d_c[0][0].size() == 4);

    // map of equation-moment order for diffusion coefficients
    d_d[0] = 1;
    d_d[1] = 3;
    d_d[2] = 5;
    d_d[3] = 7;

    // diffusion coefficients
    d_alpha[0] = fractions::rat_1_3;
    d_alpha[1] = fractions::rat_1_7;
    d_alpha[2] = fractions::rat_1_11;
    d_alpha[3] = fractions::rat_1_15;

    // map of equation-moment order for A-matrix coefficients
    d_a[0] = 0;
    d_a[1] = 2;
    d_a[2] = 4;
    d_a[3] = 6;

    // A-matrix coefficients
    d_c[0][0][0] =  1.0;
    d_c[0][1][0] = -fractions::rat_2_3;
    d_c[0][2][0] =  fractions::rat_8_15;
    d_c[0][3][0] = -fractions::rat_16_35;

    d_c[1][0][0] = -fractions::rat_2_3;
    d_c[1][1][0] =  fractions::rat_4_9;
    d_c[1][1][1] =  fractions::rat_5_9;
    d_c[1][2][0] = -fractions::rat_16_45;
    d_c[1][2][1] = -fractions::rat_4_9;
    d_c[1][3][0] =  fractions::rat_32_105;
    d_c[1][3][1] =  fractions::rat_8_21;

    d_c[2][0][0] =  fractions::rat_8_15;
    d_c[2][1][0] = -fractions::rat_16_45;
    d_c[2][1][1] = -fractions::rat_4_9;
    d_c[2][2][0] =  fractions::rat_64_225;
    d_c[2][2][1] =  fractions::rat_16_45;
    d_c[2][2][2] =  fractions::rat_9_25;
    d_c[2][3][0] = -fractions::rat_128_525;
    d_c[2][3][1] = -fractions::rat_32_105;
    d_c[2][3][2] = -fractions::rat_54_175;

    d_c[3][0][0] = -fractions::rat_16_35;
    d_c[3][1][0] =  fractions::rat_32_105;
    d_c[3][1][1] =  fractions::rat_8_21;
    d_c[3][2][0] = -fractions::rat_128_525;
    d_c[3][2][1] = -fractions::rat_32_105;
    d_c[3][2][2] = -fractions::rat_54_175;
    d_c[3][3][0] =  fractions::rat_256_1225;
    d_c[3][3][1] =  fractions::rat_64_245;
    d_c[3][3][2] =  fractions::rat_324_1225;
    d_c[3][3][3] =  fractions::rat_13_49;

    // B-matrix coefficients
    d_b[0][0] = fractions::rat_1_2;
    d_b[0][1] = -fractions::rat_1_8;
    d_b[0][2] = fractions::rat_1_16;
    d_b[0][3] = -fractions::rat_5_128;

    d_b[1][0] = -fractions::rat_1_8;
    d_b[1][1] = fractions::rat_7_24;
    d_b[1][2] = -fractions::rat_41_384;
    d_b[1][3] = fractions::rat_1_16;

    d_b[2][0] = fractions::rat_1_16;
    d_b[2][1] = -fractions::rat_41_384;
    d_b[2][2] = fractions::rat_407_1920;
    d_b[2][3] = -fractions::rat_233_2560;

    d_b[3][0] = -fractions::rat_5_128;
    d_b[3][1] = fractions::rat_1_16;
    d_b[3][2] = -fractions::rat_233_2560;
    d_b[3][3] = fractions::rat_3023_17920;

    // F-matrix coefficients
    d_f[0][0] = 1.0;
    d_f[0][1] = -fractions::rat_2_3;
    d_f[0][2] = fractions::rat_8_15;
    d_f[0][3] = -fractions::rat_16_35;

    d_f[1][0] = -fractions::rat_2_3;
    d_f[1][1] = fractions::rat_4_9;
    d_f[1][2] = -fractions::rat_16_45;
    d_f[1][3] = fractions::rat_32_105;

    d_f[2][0] = fractions::rat_8_15;
    d_f[2][1] = -fractions::rat_16_45;
    d_f[2][2] = fractions::rat_64_225;
    d_f[2][3] = -fractions::rat_128_525;

    d_f[3][0] = -fractions::rat_16_35;
    d_f[3][1] = fractions::rat_32_105;
    d_f[3][2] = -fractions::rat_128_525;
    d_f[3][3] = fractions::rat_256_1225;

    // store the minimum Pn order of cross section data
    d_min_moments = d_mat->xs().pn_order() + 1;

    // compare across all domains
    profugus::global_min(d_min_moments);

    Ensure (d_min_moments > 0);

    // Determine of Pn correction should be performed
    if( db->isParameter("Pn_correction") )
        d_outscatter_correction = db->get<bool>("Pn_correction");

    // Build hash table for Sigma
    d_Sigma = Teuchos::rcp(new Hash_Table);
    Check (!d_Sigma.is_null());

    // set the number of moments
    int num_mom = d_dim->num_moments();

    // get the material ids from the database
    Vec_Int mats;
    d_mat->xs().get_matids(mats);
    Check (mats.size() > 0);

    // iterate through materials in the database and add them to the
    // hash-table
    for(int imom = 0; imom < num_mom; ++imom)
    {
        for (Vec_Int::const_iterator m = mats.begin(); m != mats.end(); ++m)
        {
            RCP_Serial_Matrix S( new Serial_Matrix(d_Ng,d_Ng) );
            make_Sigma(imom, *m, *S);
            d_Sigma->insert(Hash_Table::value_type(to_size_type(imom, *m), S));
        }
    }

    // complete the hash-table
    d_Sigma->complete();
    Check (d_Sigma->size() == num_mom * mats.size());
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Make the \f$\boldsymbol{\Sigma}_n\f$ cross section matrix.
 *
 * The cross section matrix is a \f$N_g\times N_g\f$ matrix and is defined
 * \f[
 * \boldsymbol{\Sigma}_n = \begin{pmatrix}
 * (\sigma^0-\sigma_n^{00}) & -\sigma_n^{01} & \dots & -\sigma_n^{0G} \\
 * &&&\\
 * -\sigma_n^{10} & (\sigma^1-\sigma_n^{11}) & \dots & -\sigma_n^{1G} \\
 * &&&\\
 * \vdots & \vdots & \ddots & \vdots \\
 * &&&\\
 * -\sigma_n^{G0} & -\sigma_n^{G1} & \dots & (\sigma^G-\sigma_n^{GG})
 * \end{pmatrix}
 * \f]
 *
 * \param n moment of cross sections in range [0,7]
 * \param matid material id
 * \param S pre-allocated \f$N_g\times N_g\f$ matrix
 *
 * \pre n < N + 1 for a given \f$SP_N\f$ approximation
 */
void Moment_Coefficients::make_Sigma(int            n,
                                     int            matid,
                                     Serial_Matrix &S)
{
    Require (!d_mat.is_null());
    Require (n >= 0 && n < d_dim->num_moments());
    Require (d_mat->xs().has(matid));
    Require (S.numRows() == S.numCols());
    Require (S.numRows() == d_Ng);
    Require (d_mat->xs().num_groups() == d_Ng);

    // assign the matrix to 0.0
    S.putScalar(0.0);

    // get constant reference to the cross sections
    const XS_t &xs = d_mat->xs();

    // get the total cross sections for this mat
    const XS_t::Vector &total = xs.vector(matid, XS_t::TOTAL);
    Check (total.length() == d_Ng);

    // put the group totals on the diagonal
    for (int g = 0; g < d_Ng; ++g)
    {
        S(g, g) = total(g);
        Check (S(g, g) >= 0.0);
    }

    // add the scattering cross sections if the moments exist in the data
    if (n < d_min_moments)
    {
        // get the cross section matrix
        const XS_t::Matrix &scat = xs.matrix(matid, n);

        // loop over groups
        for (int g = 0; g < d_Ng; ++g)
        {
            for (int gp = 0; gp < d_Ng; ++gp)
            {
                S(g, gp) -= scat(g, gp);
            }
        }
    }
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Moment_Coefficients.cc
//---------------------------------------------------------------------------//
