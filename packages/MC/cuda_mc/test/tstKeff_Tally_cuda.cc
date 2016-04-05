//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstKeff_Tally_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:36:57 2016
 * \brief  Test for Keff_Tally
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>
#include <vector>

#include "Keff_Tally_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"
#include "geometry/Cartesian_Mesh.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class Keff_Tally_cudaTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::XS          XS_t;
    typedef Teuchos::RCP<XS_t>    RCP_XS;

  protected:
    void SetUp()
    {
    }

    void build_1grp_xs()
    {
        constexpr int ng = 1;
        xs = Teuchos::rcp(new XS_t());
        xs->set(0, ng);

        // make group boundaries
        XS_t::OneDArray nbnd(ng+1, 0.0);
        nbnd[0] = 100.0;
        nbnd[1] = 0.00001;
        xs->set_bounds(nbnd);

        double t1[ng] = {1.0};

        XS_t::OneDArray tot1(std::begin(t1), std::end(t1));

        xs->add(0, XS_t::TOTAL, tot1);

        double s1[][ng] = {{0.5}};

        XS_t::TwoDArray sct1(ng, ng, 0.0);

        for (int g = 0; g < ng; ++g)
        {
            for (int gp = 0; gp < ng; ++gp)
            {
                sct1(g, gp) = s1[g][gp];
            }
        }

        xs->add(0, 0, sct1);

        double c[] = {1.0};
        double f[] = {0.2};
        double n[] = {0.4};

        XS_t::OneDArray chi(std::begin(c),std::end(c));
        XS_t::OneDArray fis(std::begin(f),std::end(f));
        XS_t::OneDArray nuf(std::begin(n),std::end(n));

        xs->add(0, XS_t::CHI,       chi);
        xs->add(0, XS_t::SIG_F,     fis);
        xs->add(0, XS_t::NU_SIG_F,  nuf);

        xs->complete();
    }

    void build_2grp_xs()
    {
        constexpr int ng = 2;
        xs = Teuchos::rcp(new XS_t());
        xs->set(0, ng);

        // make group boundaries
        XS_t::OneDArray nbnd(ng+1, 0.0);
        nbnd[0] = 100.0;
        nbnd[1] = 1.0;
        nbnd[2] = 0.00001;
        xs->set_bounds(nbnd);

        double t1[ng] = {1.0, 0.1};

        XS_t::OneDArray tot1(std::begin(t1), std::end(t1));

        xs->add(0, XS_t::TOTAL, tot1);

        double s1[][ng] = {{0.5, 0.0},{0.1, 0.05}};

        XS_t::TwoDArray sct1(ng, ng, 0.0);

        for (int g = 0; g < ng; ++g)
        {
            for (int gp = 0; gp < ng; ++gp)
            {
                sct1(g, gp) = s1[g][gp];
            }
        }

        xs->add(0, 0, sct1);

        double c[] = {0.9, 0.10};
        double f[] = {0.2, 0.02};
        double n[] = {0.4, 0.05};

        XS_t::OneDArray chi(std::begin(c),std::end(c));
        XS_t::OneDArray fis(std::begin(f),std::end(f));
        XS_t::OneDArray nuf(std::begin(n),std::end(n));

        xs->add(0, XS_t::CHI,       chi);
        xs->add(0, XS_t::SIG_F,     fis);
        xs->add(0, XS_t::NU_SIG_F,  nuf);

        xs->complete();
    }

  protected:
    // >>> DATA
    RCP_XS xs;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Keff_Tally_cudaTest, onegroup)
{
    this->build_1grp_xs();

    // Mesh edges
    std::vector<double> x_edges = {0.0, 0.20, 1.0};
    std::vector<double> y_edges = {0.0, 0.50, 1.0};
    std::vector<double> z_edges = {0.0, 0.70, 1.0};

    // Cartesian mesh object to compute volumes
    profugus::Cartesian_Mesh mesh(x_edges,y_edges,z_edges);

    int num_particles = 64;
    double keff;
    cuda_mc::Keff_Tally_Tester::test_tally(x_edges,y_edges,z_edges,xs,
                                           keff,num_particles);

    std::cout << "Computed keff: " << keff << std::endl;

    // Value should be 0.6
    // nu_sigma_f is 0.4 and each thread adds contribution of 1.5
    double expected_keff = 1.5 * 0.4;
    double tol = 1.0 / std::sqrt( static_cast<double>(num_particles) );
    EXPECT_SOFTEQ( expected_keff, keff, tol );

}

TEST_F(Keff_Tally_cudaTest, twogroup)
{
    build_2grp_xs();

    // Mesh edges
    std::vector<double> x_edges = {0.0, 0.20, 1.0};
    std::vector<double> y_edges = {0.0, 0.50, 1.0};
    std::vector<double> z_edges = {0.0, 0.70, 1.0};

    // Cartesian mesh object to compute volumes
    profugus::Cartesian_Mesh mesh(x_edges,y_edges,z_edges);

    int num_particles = 64;
    double keff;
    cuda_mc::Keff_Tally_Tester::test_tally(x_edges,y_edges,z_edges,xs,
                                           keff,num_particles);

    std::cout << "Computed keff: " << keff << std::endl;

    // Value should be 0.6
    // Half threads in group 0 add 1.5 * 0.4
    // Half threads in group 1 add 1.5 * 0.05
    double expected_keff = 0.5 * (1.5 * 0.4 + 1.5 * 0.05);
    double tol = 1.0 / std::sqrt( static_cast<double>(num_particles) );
    EXPECT_SOFTEQ( expected_keff, keff, tol );

}

//---------------------------------------------------------------------------//
//                 end of tstKeff_Tally_cuda.cc
//---------------------------------------------------------------------------//
