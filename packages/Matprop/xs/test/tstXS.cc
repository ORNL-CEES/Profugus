//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/test/tstXS.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 31 12:54:58 2014
 * \brief  XS container test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../XS.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
class XS_Test : public testing::Test
{
  protected:
    typedef profugus::XS  XS;
    typedef XS::Vector    Vector;
    typedef XS::Matrix    Matrix;
    typedef XS::OneDArray OneDArray;
    typedef XS::TwoDArray TwoDArray;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        Ng = 4;

        // build data
        m1_sig.resize(Ng);

        m5_sig.resize(Ng);
        sigf.resize(Ng);
        nusigf.resize(Ng);
        chi.resize(Ng);

        m1_sig[0] = 2.0;
        m1_sig[1] = 3.0;
        m1_sig[2] = 4.0;
        m1_sig[3] = 5.0;

        m5_sig[0] = 20.0;
        m5_sig[1] = 30.0;
        m5_sig[2] = 40.0;
        m5_sig[3] = 50.0;

        sigf[0] = 11.0;
        sigf[1] = 12.0;
        sigf[2] = 13.0;
        sigf[3] = 14.0;

        double nu = 2.4;
        for (int n = 0; n < 4; ++n)
        {
            nusigf[n] = sigf[n] * nu;
        }

        chi[0] = 0.6;
        chi[1] = 0.3;
        chi[2] = 0.1;

        xs.set(1, Ng);
    }

  protected:

    int Ng;

    OneDArray m1_sig;
    OneDArray m5_sig, sigf, nusigf, chi;

    XS xs;

};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(XS_Test, totals_assignment)
{
    EXPECT_EQ(1, xs.pn_order());
    EXPECT_EQ(4, xs.num_groups());
    EXPECT_EQ(0, xs.num_mat());

    xs.add(1, XS::TOTAL, m1_sig);

    xs.add(5, XS::TOTAL, m5_sig);
    xs.add(5, XS::SIG_F, sigf);
    xs.add(5, XS::NU_SIG_F, nusigf);
    xs.add(5, XS::CHI, chi);

    EXPECT_EQ(0, xs.num_mat());

    xs.complete();
}

//---------------------------------------------------------------------------//
//                 end of tstXS.cc
//---------------------------------------------------------------------------//
