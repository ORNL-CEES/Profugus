//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/tstMoment_Coefficients.cc
 * \author Thomas M. Evans
 * \date   Mon Feb 10 14:45:24 2014
 * \brief  Test of Moment_Coefficients class.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "Teuchos_ParameterList.hpp"

#include "../Dimensions.hh"
#include "../Moment_Coefficients.hh"
#include "Test_XS.hh"

using profugus::Moment_Coefficients;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class Moment_CoefficientsTest : public testing::Test
{
  protected:
    typedef Moment_Coefficients::XS_t              XS;
    typedef Moment_Coefficients::RCP_XS            RCP_XS;
    typedef Moment_Coefficients::Mat_DB_t          Mat_DB;
    typedef Moment_Coefficients::RCP_Mat_DB        RCP_Mat_DB;
    typedef Moment_Coefficients::RCP_Dimensions    RCP_Dimensions;
    typedef Moment_Coefficients::Serial_Matrix     Serial_Matrix;
    typedef Moment_Coefficients::RCP_ParameterList RCP_ParameterList;

  protected:

    void SetUp()
    {
        mat3  = three_grp::make_mat(3, 4);
        mat12 = twelve_grp::make_mat(1, 4);

        EXPECT_FALSE(mat3.is_null());
        EXPECT_FALSE(mat12.is_null());

        db = Teuchos::rcp(new Teuchos::ParameterList("MC"));
    }

    void check_matrices(const Serial_Matrix &A,
                        const Serial_Matrix &B,
                        double eps = 1.0e-6)
    {
        int N = A.numRows();
        EXPECT_EQ(N, A.numCols());
        EXPECT_EQ(N, B.numCols());
        EXPECT_EQ(N, B.numRows());

        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < N; ++i)
            {
                EXPECT_NEAR(B(i, j), A(i, j), eps);
            }
        }
    }

    void make_dim(int N)
    {
        dim = Teuchos::rcp(new profugus::Dimensions(N));
    }

  protected:

    RCP_ParameterList db;
    RCP_Dimensions dim;
    RCP_Mat_DB mat3, mat12;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Moment_CoefficientsTest, Sigma_3grp)
{
    make_dim(7);

    Moment_Coefficients mc(db, dim, mat3);
    EXPECT_EQ(3, mc.num_groups());
    EXPECT_EQ(4, mc.min_scattering_moments());
    EXPECT_EQ(4, mc.num_equations());

    Serial_Matrix S0(3, 3), S1(3, 3), S2(3, 3), S3(3, 3), S4(3, 3),
        S5(3, 3), S6(3, 3), S7(3, 3);

    mc.make_Sigma(0, 0, S0);
    mc.make_Sigma(1, 0, S1);
    mc.make_Sigma(2, 0, S2);
    mc.make_Sigma(3, 0, S3);
    mc.make_Sigma(4, 0, S4);
    mc.make_Sigma(5, 0, S5);
    mc.make_Sigma(6, 0, S6);
    mc.make_Sigma(7, 0, S7);

    double eps = 1.0e-6;

    double r0 = 0.0, r1 = 0.0, r2 = 0.0, r3 = 0.0, r4 = 0.0, r5 = 0.0, r6 = 0.0,
           r7 = 0.0;

    // check matrices
    for (int g = 0; g < 3; ++g)
    {
        for (int gp = 0; gp < 3; ++gp)
        {
            r0 = -three_grp::S0[g][gp];
            r1 = -three_grp::S1[g][gp];
            r2 = -three_grp::S2[g][gp];
            r3 = -three_grp::S3[g][gp];
            r4 = 0.0;
            r5 = 0.0;
            r6 = 0.0;
            r7 = 0.0;

            if (g == gp)
            {
                r0 += three_grp::T[g];
                r1 += three_grp::T[g];
                r2 += three_grp::T[g];
                r3 += three_grp::T[g];
                r4 += three_grp::T[g];
                r5 += three_grp::T[g];
                r6 += three_grp::T[g];
                r7 += three_grp::T[g];
            }

            EXPECT_NEAR(r0, S0(g, gp), eps);
            EXPECT_NEAR(r1, S1(g, gp), eps);
            EXPECT_NEAR(r2, S2(g, gp), eps);
            EXPECT_NEAR(r3, S3(g, gp), eps);
            EXPECT_NEAR(r4, S4(g, gp), eps);
            EXPECT_NEAR(r5, S5(g, gp), eps);
            EXPECT_NEAR(r6, S6(g, gp), eps);
            EXPECT_NEAR(r7, S7(g, gp), eps);
        }
    }
}

//---------------------------------------------------------------------------//

TEST_F(Moment_CoefficientsTest, Sigma_12grp)
{
    make_dim(3);

    Moment_Coefficients mc(db, dim, mat12);
    EXPECT_EQ(12, mc.num_groups());
    EXPECT_EQ(2, mc.min_scattering_moments());
    EXPECT_EQ(4, dim->num_moments());

    Serial_Matrix S0(12, 12), S1(12, 12), S2(12, 12), S3(12, 12);

    mc.make_Sigma(0, 0, S0);
    mc.make_Sigma(1, 0, S1);
    mc.make_Sigma(2, 0, S2);
    mc.make_Sigma(3, 0, S3);

    double eps = 1.0e-6;

    double r0 = 0.0, r1 = 0.0, r2 = 0.0, r3 = 0.0;

    // check matrices
    for (int g = 0; g < 12; ++g)
    {
        for (int gp = 0; gp < 12; ++gp)
        {
            r0 = -twelve_grp::S0[g][gp];
            r1 = -twelve_grp::S1[g][gp];
            r2 = 0.0;
            r3 = 0.0;

            if (g == gp)
            {
                r0 += twelve_grp::T[g];
                r1 += twelve_grp::T[g];
                r2 += twelve_grp::T[g];
                r3 += twelve_grp::T[g];
            }

            EXPECT_NEAR(r0, S0(g, gp), eps);
            EXPECT_NEAR(r1, S1(g, gp), eps);
            EXPECT_NEAR(r2, S2(g, gp), eps);
            EXPECT_NEAR(r3, S3(g, gp), eps);
        }
    }
}

//---------------------------------------------------------------------------//

TEST_F(Moment_CoefficientsTest, Diffusion_3grp_SP7)
{
    make_dim(7);

    Moment_Coefficients mc(db, dim, mat3);
    EXPECT_EQ(3, mc.num_groups());
    EXPECT_EQ(4, mc.min_scattering_moments());
    EXPECT_EQ(4, dim->num_equations());

    Serial_Matrix D0(3, 3), D1(3, 3), D2(3, 3), D3(3, 3);

    mc.make_D(0, 0, D0);
    mc.make_D(1, 0, D1);
    mc.make_D(2, 0, D2);
    mc.make_D(3, 0, D3);

    double eps = 1.0e-6;

    EXPECT_NEAR(1.202353727657262, D0(0, 0), eps);
    EXPECT_NEAR(0, D0(0, 1), eps);
    EXPECT_NEAR(0, D0(0, 2), eps);
    EXPECT_NEAR(0.003745600175524, D0(1, 0), eps);
    EXPECT_NEAR(0.555595355302987, D0(1, 1), eps);
    EXPECT_NEAR(0.000028809408643, D0(1, 2), eps);
    EXPECT_NEAR(0.000047091227445, D0(2, 0), eps);
    EXPECT_NEAR(0.006985173541755, D0(2, 1), eps);
    EXPECT_NEAR(0.221595322957063, D0(2, 2), eps);

    EXPECT_NEAR(0.326781580498721, D1(0, 0), eps);
    EXPECT_NEAR(0, D1(0, 1), eps);
    EXPECT_NEAR(0, D1(0, 2), eps);
    EXPECT_NEAR(-0.000474508568926, D1(1, 0), eps);
    EXPECT_NEAR(0.138754557876912, D1(1, 1), eps);
    EXPECT_NEAR(0.000002801226602, D1(1, 2), eps);
    EXPECT_NEAR(0.000002751725326, D1(2, 0), eps);
    EXPECT_NEAR(-0.000804652341348, D1(2, 1), eps);
    EXPECT_NEAR(0.080112513018495, D1(2, 2), eps);

    EXPECT_NEAR(1.0/three_grp::T[0]/ 11.0, D2(0, 0), eps);
    EXPECT_NEAR(0, D2(0, 1), eps);
    EXPECT_NEAR(0, D2(0, 2), eps);
    EXPECT_NEAR(1.0/three_grp::T[1]/ 11.0, D2(1, 1), eps);
    EXPECT_NEAR(0, D2(1, 0), eps);
    EXPECT_NEAR(0, D2(1, 2), eps);
    EXPECT_NEAR(1.0/three_grp::T[2]/ 11.0, D2(2, 2), eps);
    EXPECT_NEAR(0, D2(2, 0), eps);
    EXPECT_NEAR(0, D2(2, 1), eps);

    EXPECT_NEAR(1.0/three_grp::T[0]/ 15.0, D3(0, 0), eps);
    EXPECT_NEAR(0, D3(0, 1), eps);
    EXPECT_NEAR(0, D3(0, 2), eps);
    EXPECT_NEAR(1.0/three_grp::T[1]/ 15.0, D3(1, 1), eps);
    EXPECT_NEAR(0, D3(1, 0), eps);
    EXPECT_NEAR(0, D3(1, 2), eps);
    EXPECT_NEAR(1.0/three_grp::T[2]/ 15.0, D3(2, 2), eps);
    EXPECT_NEAR(0, D3(2, 0), eps);
    EXPECT_NEAR(0, D3(2, 1), eps);
}

//---------------------------------------------------------------------------//

TEST_F(Moment_CoefficientsTest, Diffusion_3grp_SP5)
{
    make_dim(5);

    Moment_Coefficients mc(db, dim, mat3);
    EXPECT_EQ(3, mc.num_groups());
    EXPECT_EQ(4, mc.min_scattering_moments());
    EXPECT_EQ(3, dim->num_equations());

    Serial_Matrix D0(3, 3), D1(3, 3), D2(3, 3), D3(3, 3);

    mc.make_D(0, 0, D0);
    mc.make_D(1, 0, D1);
    mc.make_D(2, 0, D2);

    // cannot make D3 because we only have 3 equations (up to D2)
#ifdef REQUIRE_ON
    bool caught = false;
    try
    {
        mc.make_D(3, 0, D3);
    }
    catch (const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);
#endif

    double eps = 1.0e-6;

    EXPECT_NEAR(1.202353727657262, D0(0, 0), eps);
    EXPECT_NEAR(0, D0(0, 1), eps);
    EXPECT_NEAR(0, D0(0, 2), eps);
    EXPECT_NEAR(0.003745600175524, D0(1, 0), eps);
    EXPECT_NEAR(0.555595355302987, D0(1, 1), eps);
    EXPECT_NEAR(0.000028809408643, D0(1, 2), eps);
    EXPECT_NEAR(0.000047091227445, D0(2, 0), eps);
    EXPECT_NEAR(0.006985173541755, D0(2, 1), eps);
    EXPECT_NEAR(0.221595322957063, D0(2, 2), eps);

    EXPECT_NEAR(0.326781580498721, D1(0, 0), eps);
    EXPECT_NEAR(0, D1(0, 1), eps);
    EXPECT_NEAR(0, D1(0, 2), eps);
    EXPECT_NEAR(-0.000474508568926, D1(1, 0), eps);
    EXPECT_NEAR(0.138754557876912, D1(1, 1), eps);
    EXPECT_NEAR(0.000002801226602, D1(1, 2), eps);
    EXPECT_NEAR(0.000002751725326, D1(2, 0), eps);
    EXPECT_NEAR(-0.000804652341348, D1(2, 1), eps);
    EXPECT_NEAR(0.080112513018495, D1(2, 2), eps);

    EXPECT_NEAR(1.0/three_grp::T[0]/ 11.0, D2(0, 0), eps);
    EXPECT_NEAR(0, D2(0, 1), eps);
    EXPECT_NEAR(0, D2(0, 2), eps);
    EXPECT_NEAR(1.0/three_grp::T[1]/ 11.0, D2(1, 1), eps);
    EXPECT_NEAR(0, D2(1, 0), eps);
    EXPECT_NEAR(0, D2(1, 2), eps);
    EXPECT_NEAR(1.0/three_grp::T[2]/ 11.0, D2(2, 2), eps);
    EXPECT_NEAR(0, D2(2, 0), eps);
    EXPECT_NEAR(0, D2(2, 1), eps);
}

//---------------------------------------------------------------------------//

TEST_F(Moment_CoefficientsTest, A_3gp_SP7)
{
    make_dim(7);

    Moment_Coefficients mc(db, dim, mat3);
    EXPECT_EQ(3, mc.num_groups());
    EXPECT_EQ(4, mc.min_scattering_moments());
    EXPECT_EQ(4, dim->num_equations());

    Serial_Matrix A(3, 3), R(3, 3);
    Serial_Matrix S0(3, 3), S2(3, 3), S4(3, 3), S6(3, 3);

    // we have test make_Sigma, so use it here to make reference matrices

    // Row 0
    {
        mc.make_Sigma(0, 0, S0);

        mc.make_A(0, 0, 0, A);
        check_matrices(A, S0);

        mc.make_Sigma(0, 0, S0);
        S0 *= -2.0/3.0;

        mc.make_A(0, 1, 0, A);
        check_matrices(A, S0);

        mc.make_Sigma(0, 0, S0);
        S0 *= 8.0/15.0;

        mc.make_A(0, 2, 0, A);
        check_matrices(A, S0);

        mc.make_Sigma(0, 0, S0);
        S0 *= -16.0/35.0;

        mc.make_A(0, 3, 0, A);
        check_matrices(A, S0);
    }

    // Row 1
    {
        mc.make_Sigma(0, 0, S0);
        S0 *=-2.0/3.0;

        mc.make_A(1, 0, 0, A);
        check_matrices(A, S0);

        mc.make_Sigma(0, 0, S0);
        mc.make_Sigma(2, 0, S2);
        S0 *= 4.0/9.0;
        S2 *= 5.0/9.0;
        R   = S0;
        R  += S2;

        mc.make_A(1, 1, 0, A);
        check_matrices(A, R);

        mc.make_Sigma(0, 0, S0);
        mc.make_Sigma(2, 0, S2);
        S0 *= -16.0/45.0;
        S2 *= -4.0/9.0;
        R   = S0;
        R  += S2;

        mc.make_A(1, 2, 0, A);
        check_matrices(A, R);

        mc.make_Sigma(0, 0, S0);
        mc.make_Sigma(2, 0, S2);
        S0 *= 32.0/105.0;
        S2 *= 8.0/21.0;
        R   = S0;
        R  += S2;

        mc.make_A(1, 3, 0, A);
        check_matrices(A, R);
    }

    // Row 2
    {
        mc.make_Sigma(0, 0, S0);
        S0 *= 8.0/15.0;

        mc.make_A(2, 0, 0, A);
        check_matrices(A, S0);

        mc.make_Sigma(0, 0, S0);
        mc.make_Sigma(2, 0, S2);
        S0 *= -16.0/45.0;
        S2 *= -4.0/9.0;
        R   = S0;
        R  += S2;

        mc.make_A(2, 1, 0, A);
        check_matrices(A, R);

        mc.make_Sigma(0, 0, S0);
        mc.make_Sigma(2, 0, S2);
        mc.make_Sigma(4, 0, S4);
        S0 *= 64.0/225.0;
        S2 *= 16.0/45.0;
        S4 *= 9.0/25.0;
        R   = S0;
        R  += S2;
        R  += S4;

        mc.make_A(2, 2, 0, A);
        check_matrices(A, R);

        mc.make_Sigma(0, 0, S0);
        mc.make_Sigma(2, 0, S2);
        mc.make_Sigma(4, 0, S4);
        S0 *= -128.0/525.0;
        S2 *= -32.0/105.0;
        S4 *= -54.0/175.0;
        R   = S0;
        R  += S2;
        R  += S4;

        mc.make_A(2, 3, 0, A);
        check_matrices(A, R);
    }

    // Row 3
    {
        mc.make_Sigma(0, 0, S0);
        S0 *= -16.0/35.0;

        mc.make_A(3, 0, 0, A);
        check_matrices(A, S0);

        mc.make_Sigma(0, 0, S0);
        mc.make_Sigma(2, 0, S2);
        S0 *= 32.0/105.0;
        S2 *= 8.0/21.0;
        R   = S0;
        R  += S2;

        mc.make_A(3, 1, 0, A);
        check_matrices(A, R);

        mc.make_Sigma(0, 0, S0);
        mc.make_Sigma(2, 0, S2);
        mc.make_Sigma(4, 0, S4);
        S0 *= -128.0/525.0;
        S2 *= -32.0/105.0;
        S4 *= -54.0/175.0;
        R   = S0;
        R  += S2;
        R  += S4;

        mc.make_A(3, 2, 0, A);
        check_matrices(A, R);

        mc.make_Sigma(0, 0, S0);
        mc.make_Sigma(2, 0, S2);
        mc.make_Sigma(4, 0, S4);
        mc.make_Sigma(6, 0, S6);
        S0 *= 256.0/1225.0;
        S2 *= 64.0/245.0;
        S4 *= 324.0/1225.0;
        S6 *= 13.0/49.0;
        R   = S0;
        R  += S2;
        R  += S4;
        R  += S6;

        mc.make_A(3, 3, 0, A);
        check_matrices(A, R);
    }
}

//---------------------------------------------------------------------------//

TEST_F(Moment_CoefficientsTest, A_3gp_SP3)
{
    make_dim(3);

    Moment_Coefficients mc(db, dim, mat3);
    EXPECT_EQ(3, mc.num_groups());
    EXPECT_EQ(4, mc.min_scattering_moments());
    EXPECT_EQ(2, dim->num_equations());

    Serial_Matrix A(3, 3);

    // cannot make A3 because we only have 3 equations (up to D2)
#ifdef REQUIRE_ON
    bool caught = false;
    try
    {
        mc.make_A(0, 2, 0, A);
    }
    catch (const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);

    caught = false;
    try
    {
        mc.make_A(0, 3, 0, A);
    }
    catch (const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);

    caught = false;
    try
    {
        mc.make_A(2, 0, 0, A);
    }
    catch (const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);

    caught = false;
    try
    {
        mc.make_A(3, 0, 0, A);
    }
    catch (const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);
#endif
}

//---------------------------------------------------------------------------//

TEST_F(Moment_CoefficientsTest, B_3gp_SP7)
{
    make_dim(7);

    Moment_Coefficients mc(db, dim, mat3);
    EXPECT_EQ(3, mc.num_groups());
    EXPECT_EQ(4, mc.min_scattering_moments());
    EXPECT_EQ(4, dim->num_equations());

    Serial_Matrix B(3, 3), I(3, 3), R(3, 3);

    // make identity
    I.putScalar(0.0);
    I(0, 0) = 1.0;
    I(1, 1) = 1.0;
    I(2, 2) = 1.0;

    // Row 0
    {
        mc.make_B(0, 0, B);
        R  = I;
        R *= 1.0/2.0;
        check_matrices(B, R);

        mc.make_B(0, 1, B);
        R  = I;
        R *= -1.0/8.0;
        check_matrices(B, R);

        mc.make_B(0, 2, B);
        R  = I;
        R *= 1.0/16.0;
        check_matrices(B, R);

        mc.make_B(0, 3, B);
        R  = I;
        R *= -5.0/128.0;
        check_matrices(B, R);
    }

    // Row 1
    {
        mc.make_B(1, 0, B);
        R  = I;
        R *= -1.0/8.0;
        check_matrices(B, R);

        mc.make_B(1, 1, B);
        R  = I;
        R *= 7.0/24.0;
        check_matrices(B, R);

        mc.make_B(1, 2, B);
        R  = I;
        R *= -41.0/384.0;
        check_matrices(B, R);

        mc.make_B(1, 3, B);
        R  = I;
        R *= 1.0/16.0;
        check_matrices(B, R);
    }

    // Row 2
    {
        mc.make_B(2, 0, B);
        R  = I;
        R *= 1.0/16.0;
        check_matrices(B, R);

        mc.make_B(2, 1, B);
        R  = I;
        R *= -41.0/384.0;
        check_matrices(B, R);

        mc.make_B(2, 2, B);
        R  = I;
        R *= 407.0/1920.0;
        check_matrices(B, R);

        mc.make_B(2, 3, B);
        R  = I;
        R *= -233.0/2560.0;
        check_matrices(B, R);
    }

    // Row 3
    {
        mc.make_B(3, 0, B);
        R  = I;
        R *= -5.0/128.0;
        check_matrices(B, R);

        mc.make_B(3, 1, B);
        R  = I;
        R *= 1.0/16.0;
        check_matrices(B, R);

        mc.make_B(3, 2, B);
        R  = I;
        R *= -233.0/2560.0;
        check_matrices(B, R);

        mc.make_B(3, 3, B);
        R  = I;
        R *= 3023.0/17920.0;
        check_matrices(B, R);
    }
}

//---------------------------------------------------------------------------//

TEST_F(Moment_CoefficientsTest, B_3gp_SP3)
{
    make_dim(3);

    Moment_Coefficients mc(db, dim, mat3);
    EXPECT_EQ(3, mc.num_groups());
    EXPECT_EQ(4, mc.min_scattering_moments());
    EXPECT_EQ(2, dim->num_equations());

    Serial_Matrix B(3, 3);

    // cannot make A3 because we only have 3 equations (up to D2)
#ifdef REQUIRE_ON
    bool caught = false;
    try
    {
        mc.make_B(0, 2, B);
    }
    catch (const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);

    caught = false;
    try
    {
        mc.make_B(0, 3, B);
    }
    catch (const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);

    caught = false;
    try
    {
        mc.make_B(2, 0, B);
    }
    catch (const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);

    caught = false;
    try
    {
        mc.make_B(3, 0, B);
    }
    catch (const profugus::assertion &a)
    {
        caught = true;
    }
    EXPECT_TRUE(caught);
#endif
}

//---------------------------------------------------------------------------//

TEST_F(Moment_CoefficientsTest, F_3gp)
{
    make_dim(7);

    Serial_Matrix F(3, 3), R(3, 3), C(3, 3);
    R.putScalar(0.0);

    // no fission test
    {
        Moment_Coefficients mc(db, dim, mat3);

        // check with no fission
        for (int n = 0; n < 4; ++n)
        {
            for (int m = 0; m < 4; ++m)
            {
                mc.make_F(n, m, 0, F);
                check_matrices(F, R);
            }
        }
    }

    // fission test, we don't need to assign the scattering matrices because
    // they do not affect fission matrix construction
    {
        RCP_Mat_DB matf = Teuchos::rcp(new Mat_DB);
        {
            const XS &old = mat3->xs();
            RCP_XS xs     = Teuchos::rcp(new XS);

            xs->set(old.pn_order(), old.num_groups());

            XS::OneDArray tot(3), nusigf(3), chi(3);
            XS::TwoDArray P0(3, 3), P1(3, 3), P2(3, 3), P3(3, 3);

            const XS::Vector &sig = old.vector(0, XS::TOTAL);
            for (int g = 0; g < 3; ++g)
            {
                tot[g] = sig[g];
            }
            nusigf[0] = 1.3;
            nusigf[1] = 4.2;
            nusigf[2] = 0.0;

            chi[0] = 0.2;
            chi[1] = 0.8;
            chi[2] = 0.0;

            xs->add(0, XS::TOTAL, tot);
            xs->add(0, XS::NU_SIG_F, nusigf);
            xs->add(0, XS::CHI, chi);

            xs->complete();
            matf->set(xs, 4);
            for (int c = 0; c < 4; ++c)
            {
                matf->matid(c) = 0;
            }
        }

        R(0, 0) = 0.2 * 1.3;
        R(0, 1) = 0.2 * 4.2;

        R(1, 0) = 0.8 * 1.3;
        R(1, 1) = 0.8 * 4.2;

        Moment_Coefficients mc(db, dim, matf);

        // test
        mc.make_F(0, 0, 0, F);
        check_matrices(F, R);

        C = R;
        C *= 4.0/9.0;
        mc.make_F(1, 1, 0, F);
        check_matrices(F, C);

        C = R;
        C *= 64.0/225.0;
        mc.make_F(2, 2, 0, F);
        check_matrices(F, C);

        C = R;
        C *= 256.0/1225.0;
        mc.make_F(3, 3, 0, F);
        check_matrices(F, C);

        C = R;
        C *= -2.0/3.0;
        mc.make_F(0, 1, 0, F);
        check_matrices(F, C);
        mc.make_F(1, 0, 0, F);
        check_matrices(F, C);

        C = R;
        C *= 8.0/15.0;
        mc.make_F(0, 2, 0, F);
        check_matrices(F, C);
        mc.make_F(2, 0, 0, F);
        check_matrices(F, C);

        C = R;
        C *= -16.0/35.0;
        mc.make_F(0, 3, 0, F);
        check_matrices(F, C);
        mc.make_F(3, 0, 0, F);
        check_matrices(F, C);

        C = R;
        C *= -16.0/45.0;
        mc.make_F(1, 2, 0, F);
        check_matrices(F, C);
        mc.make_F(2, 1, 0, F);
        check_matrices(F, C);

        C = R;
        C *= 32.0/105.0;
        mc.make_F(1, 3, 0, F);
        check_matrices(F, C);
        mc.make_F(3, 1, 0, F);
        check_matrices(F, C);

        C = R;
        C *= -128.0/525.0;
        mc.make_F(2, 3, 0, F);
        check_matrices(F, C);
        mc.make_F(3, 2, 0, F);
        check_matrices(F, C);
    }
}

//---------------------------------------------------------------------------//

TEST(Static_Functions, Convert_U_to_Phi)
{
    double u0 = 11.3;
    double u1 = -0.1;
    double u2 = 4.2;
    double u3 = 1.1;

    double phi = Moment_Coefficients::u_to_phi(u0, u1, u2, u3);

    EXPECT_SOFTEQ(13.103809523809524, phi, 1.0e-12);
}

//---------------------------------------------------------------------------//
//                 end of tstMoment_Coefficients.cc
//---------------------------------------------------------------------------//
