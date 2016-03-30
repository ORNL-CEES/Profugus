//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstPhysics_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:36:57 2016
 * \brief  Test for Physics
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>
#include <vector>

#include "Physics_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class Physics_cudaTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::XS          XS_t;
    typedef Teuchos::RCP<XS_t>    RCP_XS;

  protected:
    void SetUp()
    {
        build_xs();
    }

    void build_xs()
    {
        xs = Teuchos::rcp(new XS_t());
        xs->set(0, 5);

        std::vector<double> bnd(6, 0.0);
        bnd[0] = 100.0;
        bnd[1] = 10.0;
        bnd[2] = 1.0;
        bnd[3] = 0.1;
        bnd[4] = 0.01;
        bnd[5] = 0.001;
        xs->set_bounds(bnd);

        typename XS_t::OneDArray total(5);
        typename XS_t::TwoDArray scat(5, 5);

        double f[5] = {0.1, 0.4, 1.8, 5.7, 9.8};
        double c[5] = {0.3770, 0.4421, 0.1809, 0.0, 0.0};
        double n[5] = {2.4*f[0], 2.4*f[1], 2.4*f[2], 2.4*f[3], 2.4*f[4]};
        typename XS_t::OneDArray fission(std::begin(f), std::end(f));
        typename XS_t::OneDArray chi(std::begin(c), std::end(c));
        typename XS_t::OneDArray nus(std::begin(n), std::end(n));
        xs->add(1, XS_t::SIG_F, fission);
        xs->add(1, XS_t::NU_SIG_F, nus);
        xs->add(1, XS_t::CHI, chi);

        // mat 0
        total[0] = 5.2 ;
        total[1] = 11.4;
        total[2] = 18.2;
        total[3] = 29.9;
        total[4] = 27.3;
        xs->add(0, XS_t::TOTAL, total);

        // mat 1
        total[0] = 5.2  + f[0];
        total[1] = 11.4 + f[1];
        total[2] = 18.2 + f[2];
        total[3] = 29.9 + f[3];
        total[4] = 27.3 + f[4];
        xs->add(1, XS_t::TOTAL, total);

        scat(0, 0) = 1.2;
        scat(1, 0) = 0.9;
        scat(1, 1) = 3.2;
        scat(2, 0) = 0.4;
        scat(2, 1) = 2.8;
        scat(2, 2) = 6.9;
        scat(2, 3) = 1.5;
        scat(3, 0) = 0.1;
        scat(3, 1) = 2.1;
        scat(3, 2) = 5.5;
        scat(3, 3) = 9.7;
        scat(3, 4) = 2.1;
        scat(4, 1) = 0.2;
        scat(4, 2) = 1.3;
        scat(4, 3) = 6.6;
        scat(4, 4) = 9.9;
        xs->add(0, 0, scat);
        xs->add(1, 0, scat);

        xs->complete();
    }

  protected:
    // >>> DATA
    RCP_XS xs;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Physics_cudaTest, total)
{
    // Mesh edges
    std::vector<double> edges = {0.0, 0.50, 1.0};
    std::vector<int> matids = {0, 1, 1, 0, 0, 1, 1, 0};

    int num_p = 16;
    std::vector<double> totals(num_p);
    cuda_mc::Physics_Tester::test_total(edges,edges,edges,matids,xs,totals);

    for( int i = 0; i < num_p; ++i )
    {
        int g = i % 5;
        int matid = i % 2;
        auto expected = xs->vector(matid,profugus::XS::TOTAL);
        EXPECT_SOFT_EQ( expected[g], totals[i] );
    }
}

TEST_F(Physics_cudaTest, collide)
{
    // Mesh edges
    std::vector<double> edges = {0.0, 0.50, 1.0};
    std::vector<int> matids = {0, 1, 1, 0, 0, 1, 1, 0};

    int num_p = 16;
    //std::vector<double> totals(num_p);
    cuda_mc::Physics_Tester::test_collide(edges,edges,edges,matids,xs,num_p);

    /*
    for( int i = 0; i < num_p; ++i )
    {
        int g = i % 5;
        int matid = i % 2;
        auto expected = xs->vector(matid,profugus::XS::TOTAL);
        EXPECT_SOFT_EQ( expected[g], totals[i] );
    }
    */
}

//---------------------------------------------------------------------------//
//                 end of tstPhysics_cuda.cc
//---------------------------------------------------------------------------//
