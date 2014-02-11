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
    typedef Moment_Coefficients::RCP_Mat_DB        RCP_Mat_DB;
    typedef Moment_Coefficients::RCP_Dimensions    RCP_Dimensions;
    typedef Moment_Coefficients::Serial_Matrix     Serial_Matrix;
    typedef Moment_Coefficients::RCP_ParameterList RCP_ParameterList;

  protected:

    void SetUp()
    {
        //mat3  = three_grp::make_mat(3, 4);
        mat12 = twelve_grp::make_mat(1, 4);

        //EXPECT_TRUE(mat3);
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
//                 end of tstMoment_Coefficients.cc
//---------------------------------------------------------------------------//
