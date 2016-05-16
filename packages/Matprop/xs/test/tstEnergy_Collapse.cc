//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Matprop/xs/test/tstEnergy_Collapse.cc
 * \author Steven Hamilton
 * \date   Thu Mar 28 13:26:47 2013
 * \brief  Energy_Collapse unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include <vector>
#include <cmath>
#include <sstream>

#include "Teuchos_RCP.hpp"

#include "utils/Definitions.hh"
#include "../Mat_DB.hh"
#include "../Energy_Collapse.hh"

//---------------------------------------------------------------------------//
// TEST FIXTURE
//---------------------------------------------------------------------------//

class Energy_Collapse_Test : public testing::Test
{
  protected:
    // Typedefs.
    typedef profugus::Energy_Collapse   Energy_Collapse;
    typedef Energy_Collapse::RCP_Mat_DB RCP_Mat_DB;
    typedef Energy_Collapse::Mat_DB_t   Mat_DB_t;
    typedef Mat_DB_t::XS_t              XS_t;
    typedef Mat_DB_t::RCP_XS            RCP_XS;
    typedef def::Vec_Dbl                Vec_Dbl;
    typedef Mat_DB_t::Vec_Int           Vec_Int;
    typedef XS_t::OneDArray             OneDArray;
    typedef XS_t::TwoDArray             TwoDArray;
    typedef XS_t::Vector                Vector;
    typedef XS_t::Matrix                Matrix;

  protected:
    void SetUp()
    {
        // make cross sections
        xs = Teuchos::rcp(new XS_t);
        xs->set(0, 8);

        // 2 materials
        OneDArray tot0(8, 0.0), tot1(8, 0.0);
        for (int g = 0; g < 8; ++g)
        {
            tot0[g] = 10.0*g + 1.0;
            tot1[g] = 20.0*g + 1.0;
        }

        TwoDArray sct0(8, 8, 0.0), sct1(8, 8, 0.0);
        for (int g = 0; g < 8; ++g)
        {
            for (int gp = 0; gp <=g; ++gp)
            {
                sct0(g, gp) = static_cast<double>(g) + gp*0.01;
                sct1(g, gp) = static_cast<double>(g) + gp*0.001;
            }
        }

        // Now a few upscattering entries
        sct0(3, 4) = 3.04;
        sct0(5, 7) = 5.07;
        sct0(6, 7) = 6.07;

        sct1(4, 5) = 4.005;
        sct1(5, 6) = 5.006;
        sct1(5, 7) = 5.007;
        sct1(6, 7) = 6.007;

        xs->add(0, XS_t::TOTAL, tot0);
        xs->add(1, XS_t::TOTAL, tot1);

        xs->add(0, 0, sct0);
        xs->add(1, 0, sct1);

        xs->complete();
    }

  protected:
    RCP_XS xs;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Energy_Collapse_Test, Multi_Material)
{
    // make a mat db (for 10 cells)
    RCP_Mat_DB mat_db = Teuchos::rcp(new Mat_DB_t);
    mat_db->set(xs, 10);

    Vec_Int matids(10, 0);
    matids[1] = 1;
    matids[2] = 1;
    matids[6] = 1;
    matids[9] = 1;

    mat_db->assign(matids);

    // Build coarse mat_db
    Vec_Int steer(3);
    steer[0] = 3;
    steer[1] = 2;
    steer[2] = 3;

    Vec_Dbl weights(8);
    weights[0] = 1.0;
    weights[1] = 2.0;
    weights[2] = 3.0;
    weights[3] = 2.0;
    weights[4] = 1.0;
    weights[5] = 1.0;
    weights[6] = 2.0;
    weights[7] = 1.0;
    RCP_Mat_DB coarse_mat = Energy_Collapse::collapse_all_mats(
        mat_db, steer, weights);

    // Make sure matids got assigned correctly
    EXPECT_EQ(0, coarse_mat->matid(0));
    EXPECT_EQ(0, coarse_mat->matid(3));
    EXPECT_EQ(0, coarse_mat->matid(4));
    EXPECT_EQ(0, coarse_mat->matid(5));
    EXPECT_EQ(0, coarse_mat->matid(7));
    EXPECT_EQ(0, coarse_mat->matid(8));

    EXPECT_EQ(1, coarse_mat->matid(1));
    EXPECT_EQ(1, coarse_mat->matid(2));
    EXPECT_EQ(1, coarse_mat->matid(6));
    EXPECT_EQ(1, coarse_mat->matid(9));

    // get cross sections
    const XS_t &xs = coarse_mat->xs();

    EXPECT_EQ(3, xs.num_groups());
    EXPECT_EQ(0, xs.pn_order());

    // check mat0
    {
        const Vector &tot = xs.vector(0, XS_t::TOTAL);
        const Matrix &sct = xs.matrix(0, 0);

        EXPECT_DOUBLE_EQ( 86.0/6.0, tot[0]);
        EXPECT_DOUBLE_EQ(103.0/3.0, tot[1]);
        EXPECT_DOUBLE_EQ(244.0/4.0, tot[2]);

        EXPECT_DOUBLE_EQ(15.1/6.0, sct(0, 0));
        EXPECT_DOUBLE_EQ(0.0,      sct(0, 1));
        EXPECT_DOUBLE_EQ(0.0,      sct(0, 2));

        EXPECT_DOUBLE_EQ(42.16/6.0, sct(1, 0));
        EXPECT_DOUBLE_EQ(21.2/3.0,  sct(1, 1));
        EXPECT_DOUBLE_EQ(0.0,       sct(1, 2));

        EXPECT_DOUBLE_EQ(108.24/6.0, sct(2, 0));
        EXPECT_DOUBLE_EQ(54.3/3.0,   sct(2, 1));
        EXPECT_DOUBLE_EQ(62.6/4.0,   sct(2, 2));
    }

    // check mat1
    {
        const Vector &tot = xs.vector(1, XS_t::TOTAL);
        const Matrix &sct = xs.matrix(1, 0);

        EXPECT_DOUBLE_EQ(166.0/6.0, tot[0]);
        EXPECT_DOUBLE_EQ(203.0/3.0, tot[1]);
        EXPECT_DOUBLE_EQ(484.0/4.0, tot[2]);

        EXPECT_DOUBLE_EQ(15.01/6.0, sct(0, 0));
        EXPECT_DOUBLE_EQ(0.0,       sct(0, 1));
        EXPECT_DOUBLE_EQ(0.0,       sct(0, 2));

        EXPECT_DOUBLE_EQ(42.016/6.0, sct(1, 0));
        EXPECT_DOUBLE_EQ(18.016/3.0, sct(1, 1));
        EXPECT_DOUBLE_EQ(4.005/4.0,  sct(1, 2));

        EXPECT_DOUBLE_EQ(108.024/6.0, sct(2, 0));
        EXPECT_DOUBLE_EQ(54.03/3.0,   sct(2, 1));
        EXPECT_DOUBLE_EQ(72.072/4.0,  sct(2, 2));
    }
}

//---------------------------------------------------------------------------//
//                        end of tstEnergy_Collapse.cc
//---------------------------------------------------------------------------//
