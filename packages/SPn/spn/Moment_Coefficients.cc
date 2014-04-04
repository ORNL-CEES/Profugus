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
                                         RCP_Mat_DB        mat,
                                         RCP_Timestep      dt)


    : d_dim(dim)
    , d_mat(mat)
    , d_dt(dt)
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
    const Vector &total = xs.vector(matid, XS_t::TOTAL);
    Check (total.length() == d_Ng);

    // put the group totals on the diagonal
    for (int g = 0; g < d_Ng; ++g)
    {
        S(g, g) = total(g);
        Check (S(g, g) >= 0.0);
    }

    // add 1/vdT to diagonal for time-dependent problems
    if (!d_dt.is_null())
    {
        // get group velocities
        const auto &v = d_mat->xs().velocities();
        Check (v.length() == d_Ng);

        double inv_dt = 1.0 / d_dt->dt();

        // put the group totals on the diagonal
        for (int g = 0; g < d_Ng; ++g)
        {
            Check (v(g) > 0.0);
            S(g, g) += 1.0 / v(g) * inv_dt;
            Check (S(g, g) >= 0.0);
        }
    }

    // add the scattering cross sections if the moments exist in the data
    if (n < d_min_moments)
    {
        // get the cross section matrix
        const Matrix &scat = xs.matrix(matid, n);

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

//---------------------------------------------------------------------------//
/*!
 * \brief Make the diffusion matrix.
 *
 * The diffusion matrix is defined
 * \f[
   \mathbf{D}_n = \alpha_n(\boldsymbol{\Sigma}_{d(n)})^{-1}
 * \f]
 * where
 * \f[
   d(0) = 1\:,\quad d(1) = 3\:,\quad d(2) = 5\:, \quad d(3) = 7
 * \f]
 * and
 * \f[
   \alpha_0 = 1/3\:,\quad
   \alpha_1 = 1/7\:,\quad
   \alpha_2 = 1/11\:,\quad
   \alpha_3 = 1/15
 * \f]
 *
 * \param n equation-order in range [0,4)
 * \param cell problem cell
 * \param D pre-allocated \f$N_g\times N_g\f$ matrix
 *
 * \pre n < (N + 1)/2 for a given \f$SP_N\f$ approximation
 */
void Moment_Coefficients::make_D(int            n,
                                 int            cell,
                                 Serial_Matrix &D)
{
    Require (!d_mat.is_null());
    Require (n >= 0 && n < d_dim->num_equations());
    Require (cell < d_mat->num_cells());
    Require (D.numRows() == D.numCols());
    Require (D.numRows() == d_Ng);
    Require (d_mat->xs().num_groups() == d_Ng);

    // first get sigma for this diffusion coefficient
    int matid = d_mat->matid(cell);
    Check( d_Sigma->exists(to_size_type(d_d[n],matid)) );
    Teuchos::RCP<Serial_Matrix> S = d_Sigma->at(to_size_type(d_d[n],matid));
    Check( !S.is_null() );
    D.assign(*S);

    if( !d_outscatter_correction )
    {
        // LU decomposition
        d_lapack.GETRF(d_Ng, d_Ng, D.values(), D.stride(), &d_ipiv[0],
                       &d_info);
        Check (d_info == 0);

        // inverse
        d_lapack.GETRI(d_Ng, D.values(), D.stride(), &d_ipiv[0], &d_work[0],
                       d_Ng, &d_info);
        Check (d_info == 0);
    }
    else
    {
        Serial_Matrix sig(D);
        D.putScalar(0.0);

        // Apply outscatter correction
        for( int ig=0; ig<d_Ng; ++ig )
        {
            for( int jg=0; jg<d_Ng; ++jg )
            {
                D(ig,ig) += sig(jg,ig);
            }
        }

        // Invert diagonal entries
        for( int ig=0; ig<d_Ng; ++ig )
        {
            D(ig,ig) = 1.0/D(ig,ig);
        }
    }

    // multiply by the scalar coefficient to complete the diffusion matrix
    // definition
    D *= d_alpha[n];
}

//---------------------------------------------------------------------------//
/*
 * \brief Make A-matrix entries.
 *
 * The A matrix is a \f$4\times 4\f$ block matrix where each block is
 * \f$N_g\times N_g\f$.  The A-matrix is defined
 * \f[
   \ve{A} = \begin{pmatrix}
    (\boldsymbol{\Sigma}_0) &
    (-\frac{2}{3}\boldsymbol{\Sigma}_0) &
    (\frac{8}{15}\boldsymbol{\Sigma}_0) &
    (-\frac{16}{35}\boldsymbol{\Sigma}_0) \\
    %%
    &&&\\
    %%
    (-\frac{2}{3}\boldsymbol{\Sigma}_0) &
    (\frac{4}{15}\boldsymbol{\Sigma}_0 + \frac{1}{3}\boldsymbol{\Sigma}_2) &
    (-\frac{16}{45}\boldsymbol{\Sigma}_0 - \frac{4}{9}\boldsymbol{\Sigma}_2) &
    (\frac{32}{105}\boldsymbol{\Sigma}_0 + \frac{8}{21}\boldsymbol{\Sigma}_2) \ \
    %%
    &&&\\
    %%
    (\frac{8}{15}\boldsymbol{\Sigma}_0) &
    (-\frac{16}{45}\boldsymbol{\Sigma}_0 - \frac{4}{9}\boldsymbol{\Sigma}_2) &
    (\frac{64}{225}\boldsymbol{\Sigma}_0 + \frac{16}{45}\boldsymbol{\Sigma}_2
    + \frac{9}{25}\boldsymbol{\Sigma}_4) &
    (-\frac{128}{525}\boldsymbol{\Sigma}_0 -
    \frac{32}{105}\boldsymbol{\Sigma}_2 - \frac{54}{175}\boldsymbol{\Sigma}_4)
    \\
    %%
    &&&\\
    %%
    (-\frac{16}{35}\boldsymbol{\Sigma}_0) &
    (\frac{32}{105}\boldsymbol{\Sigma}_0 +
    \frac{8}{21}\boldsymbol{\Sigma}_2) &
    (-\frac{128}{525}\boldsymbol{\Sigma}_0 -
    \frac{32}{105}\boldsymbol{\Sigma}_2 -
    \frac{54}{175}\boldsymbol{\Sigma}_4) &
    (\frac{256}{1225}\boldsymbol{\Sigma}_0 +
    \frac{64}{245}\boldsymbol{\Sigma}_2 +
    \frac{324}{1225}\boldsymbol{\Sigma}_4 +
    \frac{13}{49}\boldsymbol{\Sigma}_6)
  \end{pmatrix}
 * \f]
 *
 * \param n row of A-matrix in range [0,4)
 * \param m column of A-matrix in range [0,4)
 * \param cell problem cell
 * \param A pre-allocated \f$N_g\times N_g\f$ matrix
 *
 * \pre n,m < (N + 1)/2 for a given \f$SP_N\f$ approximation
 */
void Moment_Coefficients::make_A(int            n,
                                 int            m,
                                 int            cell,
                                 Serial_Matrix &A)
{
    Require (!d_mat.is_null());
    Require (n >= 0 && n < d_dim->num_equations());
    Require (m >= 0 && m < d_dim->num_equations());
    Require (cell < d_mat->num_cells());
    Require (A.numRows() == A.numCols());
    Require (A.numRows() == d_Ng);
    Require (d_mat->xs().num_groups() == d_Ng);

    // initialize A to the first term in each series entry (Sigma_0)
    int matid = d_mat->matid(cell);

    Check( d_Sigma->exists(to_size_type(0,matid)) );
    Teuchos::RCP<Serial_Matrix> S = d_Sigma->at(to_size_type(0,matid));
    Check( !S.is_null() );
    A.assign(*S);

    A *= d_c[n][m][0];

    // loop over all possible linear combinations for each element in the
    // matrix
    for (int k = 1; k < 4; ++k)
    {
        // only go to the effort for non-zero entries
        if (std::fabs(d_c[n][m][k]) > 0.0)
        {
            // make sigma for this iterate in a work matrix
            Check( d_Sigma->exists(to_size_type(d_a[k],matid)) );
            S = d_Sigma->at(to_size_type(d_a[k],matid));
            Check( !S.is_null() );
            d_W = *S;

            // multiply by the scalar coefficient
            d_W *= d_c[n][m][k];

            // add it to the running total
            A += d_W;
        }
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Make B-matrix entries.
 *
 * The B matrix is a \f$4\times 4\f$ block matrix where each block is
 * \f$N_g\times N_g\f$ with non-zero entries only on the diagonal.  The
 * B-matrix is defined:
 * \f[
   \ve{B} = \begin{pmatrix}
    \frac{1}{2} &
    -\frac{1}{8} &
    \frac{1}{16} &
    -\frac{5}{128} \\
    %%
    &&&\\
    %%
    -\frac{1}{8} &
    \frac{7}{24} &
    -\frac{41}{384} &
    \frac{1}{16} \\
    %%
    &&&\\
    %%
    \frac{1}{16} &
    -\frac{41}{384} &
    \frac{407}{1920} &
    -\frac{233}{2560} \\
    %%
    &&&\\
    %%
    -\frac{5}{128} &
    \frac{1}{16} &
    -\frac{233}{2560} &
    \frac{3023}{17920}
  \end{pmatrix}
 * \f]
 * Thus, each block in the B-matrix is a diagonal matrix with a constant value
 * on the diagonal.  Furthermore, there is no spatial dependence on the
 * B-matrix.
 *
 * \param n row of B-matrix in range [0,4)
 * \param m column of B-matrix in range [0,4)
 * \param B pre-allocated \f$N_g\times N_g\f$ matrix
 *
 * \pre n,m < (N + 1)/2 for a given \f$SP_N\f$ approximation
 */
void Moment_Coefficients::make_B(int            n,
                                 int            m,
                                 Serial_Matrix &B)
{
    Require (n >= 0 && n < d_dim->num_equations());
    Require (m >= 0 && m < d_dim->num_equations());
    Require (B.numRows() == B.numCols());
    Require (B.numRows() == d_Ng);

    // add the appropriate cofficient to the diagonal
    B.putScalar(0.0);
    for (int g = 0; g < d_Ng; ++g)
    {
        B(g, g) = d_b[n][m];
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get F fission matrix block entries.
 * \brief Make the \f$\boldsymbol{\Sigma}_n\f$ cross section matrix.
 *
 * The F matrix is a \f$4\times 4\f$ block matrix where each block is
 * \f$N_g\times N_g\f$.  The F-matrix is defined
 * \f[
  \mathbf{F} = \begin{pmatrix}
    \boldsymbol{\Sigma}_{\mathrm{f}} &
    -\frac{2}{3}\boldsymbol{\Sigma}_{\mathrm{f}} &
    \frac{8}{15}\boldsymbol{\Sigma}_{\mathrm{f}} &
    -\frac{16}{35}\boldsymbol{\Sigma}_{\mathrm{f}} \\
    &&&\\
    -\frac{2}{3}\boldsymbol{\Sigma}_{\mathrm{f}} &
    \frac{4}{9}\boldsymbol{\Sigma}_{\mathrm{f}} &
    -\frac{16}{45}\boldsymbol{\Sigma}_{\mathrm{f}} &
    \frac{32}{105}\boldsymbol{\Sigma}_{\mathrm{f}} \\
    &&&\\
    \frac{8}{15}\boldsymbol{\Sigma}_{\mathrm{f}} &
    -\frac{16}{45}\boldsymbol{\Sigma}_{\mathrm{f}} &
    \frac{64}{225}\boldsymbol{\Sigma}_{\mathrm{f}} &
    -\frac{128}{525}\boldsymbol{\Sigma}_{\mathrm{f}} \\
    &&&\\
    -\frac{16}{35}\boldsymbol{\Sigma}_{\mathrm{f}} &
    \frac{32}{105}\boldsymbol{\Sigma}_{\mathrm{f}} &
    -\frac{128}{525}\boldsymbol{\Sigma}_{\mathrm{f}} &
    \frac{256}{1225}\boldsymbol{\Sigma}_{\mathrm{f}}
  \end{pmatrix}
 * \f]
 * Where the block fission matrices are
 * \f[
  \boldsymbol{\Sigma}_{\mathrm{f}} = \begin{pmatrix}
    \chi^0\nu\sigma_f^0 & \chi^0\nu\sigma_f^1  & \dots & \chi^0\nu\sigma_f^G \\
    &&&\\
    \chi^1\nu\sigma_f^0 & \chi^1\nu\sigma_f^1  & \dots & \chi^1\nu\sigma_f^G \\
    &&&\\
    \vdots & \vdots & \ddots & \vdots \\
    &&&\\
    \chi^G\nu\sigma_f^0 & \chi^G\nu\sigma_f^1  & \dots & \chi^G\nu\sigma_f^G
  \end{pmatrix}
 * \f]
 *
 * \param n row of B-matrix in range [0,4)
 * \param m column of B-matrix in range [0,4)
 * \param cell problem cell
 * \param F pre-allocated \f$N_g\times N_g\f$ matrix
 */
void Moment_Coefficients::make_F(int            n,
                                 int            m,
                                 int            cell,
                                 Serial_Matrix &F)
{
    Require (!d_mat.is_null());
    Require (n >= 0 && n < d_dim->num_equations());
    Require (m >= 0 && m < d_dim->num_equations());
    Require (cell < d_mat->num_cells());
    Require (F.numRows() == F.numCols());
    Require (F.numRows() == d_Ng);
    Require (d_mat->xs().num_groups() == d_Ng);

    // initialize F
    F.putScalar(0.0);

    // f*chi for each group
    double fchi = 0.0;

    // matid
    int matid = d_mat->matid(cell);

    // cross sections
    const XS_t &xs = d_mat->xs();

    // get the fission data (these will be full of zeros if there is no
    // fission in this material)
    const Vector &nusigf = xs.vector(matid, XS_t::NU_SIG_F);
    const Vector &chi    = xs.vector(matid, XS_t::CHI);
    Check (nusigf.length() == d_Ng);
    Check (chi.length() == d_Ng);

    // loop over the rows of F
    for (int g = 0; g < d_Ng; ++g)
    {
        // calculate chi for this row (group)
        fchi = d_f[n][m] * chi[g];

        // loop through the columns in this row
        for (int gp = 0; gp < d_Ng; ++gp)
        {
            // add the matrix entry
            F(g, gp) = fchi * nusigf(gp);
        }
    }
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Moment_Coefficients.cc
//---------------------------------------------------------------------------//
