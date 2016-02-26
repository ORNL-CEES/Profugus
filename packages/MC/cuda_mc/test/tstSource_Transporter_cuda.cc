//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstSource_Transporter_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:36:57 2016
 * \brief  Test for Source_Transporter
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>
#include <vector>

#include "Source_Transporter_Tester.hh"
#include "mc/Definitions.hh"

#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class Source_Transporter_cudaTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::XS          XS_t;
    typedef std::shared_ptr<XS_t> SP_XS;

  protected:
    void SetUp()
    {
    }

    void build_5grp_xs()
    {
        xs = SP_XS(new XS_t());
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

    void build_3grp_xs()
    {
        const int ng = 3;
        xs = SP_XS(new XS_t());
        xs->set(0, ng);

        // make group boundaries
        XS_t::OneDArray nbnd(ng+1, 0.0);
        nbnd[0] = 100.0;
        nbnd[1] = 1.0;
        nbnd[2] = 0.01;
        nbnd[3] = 0.0001;
        xs->set_bounds(nbnd);

        double t1[ng] = {1.1, 1.6, 2.9};
        double t2[ng] = {10.0, 11.3, 16.2};

        XS_t::OneDArray tot1(std::begin(t1), std::end(t1));
        XS_t::OneDArray tot2(std::begin(t2), std::end(t2));

        xs->add(0, XS_t::TOTAL, tot1);
        xs->add(1, XS_t::TOTAL, tot2);

        double s1[][3] = {{0.7, 0.0, 0.0},
                          {0.2, 0.3, 0.0},
                          {0.1, 0.7, 1.9}};

        double s2[][3] = {{2.7, 0.0, 0.0},
                          {2.2, 2.3, 0.0},
                          {2.1, 2.7, 3.9}};

        XS_t::TwoDArray sct1(ng, ng, 0.0);
        XS_t::TwoDArray sct2(ng, ng, 0.0);

        for (int g = 0; g < 3; ++g)
        {
            for (int gp = 0; gp < 3; ++gp)
            {
                sct1(g, gp) = s1[g][gp];
                sct2(g, gp) = s2[g][gp];
            }
        }

        xs->add(0, 0, sct1);
        xs->add(1, 0, sct2);

        double c2[] = {0.4, 0.6, 0.0};
        double f2[] = {3.2, 4.2, 0.0};
        double n2[] = {2.4*3.2, 2.4*4.2, 0.0};

        XS_t::OneDArray chi2(std::begin(c2), std::end(c2));
        XS_t::OneDArray fis2(std::begin(f2), std::end(f2));
        XS_t::OneDArray nuf2(std::begin(n2), std::end(n2));

        xs->add(1, XS_t::CHI, chi2);
        xs->add(1, XS_t::SIG_F, fis2);
        xs->add(1, XS_t::NU_SIG_F, nuf2);

        xs->complete();
    }

    void build_1grp1mat_xs()
    {
        const int ng = 1;
        xs = SP_XS(new XS_t());
        xs->set(0, ng);

        // make group boundaries
        XS_t::OneDArray nbnd(ng+1, 0.0);
        nbnd[0] = 100.0;
        nbnd[1] = 0.00001;
        xs->set_bounds(nbnd);

        double t1[ng] = {1.0};

        XS_t::OneDArray tot1(std::begin(t1), std::end(t1));

        xs->add(0, XS_t::TOTAL, tot1);

        double s1[][3] = {{0.5}};

        XS_t::TwoDArray sct1(ng, ng, 0.0);

        for (int g = 0; g < ng; ++g)
        {
            for (int gp = 0; gp < ng; ++gp)
            {
                sct1(g, gp) = s1[g][gp];
            }
        }

        xs->add(0, 0, sct1);

        xs->complete();
    }

  protected:

    // >>> DATA
    SP_XS xs;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Source_Transporter_cudaTest, five_group)
{
    build_5grp_xs();

    // Mesh edges
    std::vector<double> edges = {0.0, 0.50, 1.0};
    std::vector<unsigned int> matids = {0, 1, 1, 0, 0, 1, 1, 0};

    int num_p = 128;
    cuda_mc::Source_Transporter_Tester::test_transport(
        edges,edges,edges,matids,xs,num_p);

    // No tallies yet...
}

TEST_F(Source_Transporter_cudaTest, three_group)
{
    build_3grp_xs();

    // Mesh edges
    std::vector<double> edges = {0.0, 0.5, 1.0};
    std::vector<unsigned int> matids = {0, 1, 1, 0, 0, 1, 1, 0};

    int num_p = 128;
    cuda_mc::Source_Transporter_Tester::test_transport(
        edges,edges,edges,matids,xs,num_p);

    // No tallies yet...
}

TEST_F(Source_Transporter_cudaTest, one_group)
{
    build_1grp1mat_xs();

    // Mesh edges
    std::vector<double> edges = {0.0, 5.0, 10.0};
    std::vector<unsigned int> matids = {0, 0, 0, 0, 0, 0, 0, 0};

    int num_p = 128;
    cuda_mc::Source_Transporter_Tester::test_transport(
        edges,edges,edges,matids,xs,num_p);

    // No tallies yet...
}

//---------------------------------------------------------------------------//
//                 end of tstSource_Transporter_cuda.cc
//---------------------------------------------------------------------------//
