//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstCell_Tally_cuda.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:36:57 2016
 * \brief  Test for Cell_Tally
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>
#include <vector>

#include "Cell_Tally_Tester.hh"

#include "Utils/gtest/utils_gtest.hh"
#include "geometry/Cartesian_Mesh.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class Cell_Tally_cudaTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::XS          XS_t;
    typedef std::shared_ptr<XS_t> SP_XS;

  protected:
    void SetUp()
    {
        build_xs();
    }

    void build_xs()
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

TEST_F(Cell_Tally_cudaTest, tally)
{
    // Mesh edges
    std::vector<double> x_edges = {0.0, 0.20, 1.0};
    std::vector<double> y_edges = {0.0, 0.50, 1.0};
    std::vector<double> z_edges = {0.0, 0.70, 1.0};

    // Cartesian mesh object to compute volumes
    profugus::Cartesian_Mesh mesh(x_edges,y_edges,z_edges);
    std::cout << "Cell volumes: ";
    for( int cell = 0; cell < mesh.num_cells(); ++cell )
        std::cout << mesh.volume(cell) << " ";
    std::cout << std::endl;

    int num_particles = 8192;
    std::vector<int> cells = {0, 1, 2, 4, 7};
    std::vector<double> tally;
    cuda_mc::Cell_Tally_Tester::test_tally(x_edges,y_edges,z_edges,xs,
                                           cells,tally,num_particles);

    std::cout << "Tally result: ";
    for( auto t : tally )
        std::cout << t << " ";
    std::cout << std::endl;

    EXPECT_EQ( tally.size(), cells.size() );

    // Each value should be 1.5 within statistical noise
    std::vector<double> expected(cells.size(),1.5);
    double tol = 10.0 / std::sqrt( static_cast<double>(num_particles) );
    EXPECT_VEC_SOFTEQ( expected, tally, tol );

}

//---------------------------------------------------------------------------//
//                 end of tstCell_Tally_cuda.cc
//---------------------------------------------------------------------------//
