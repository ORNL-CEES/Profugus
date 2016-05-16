//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/test/tstGeneral_Source.cc
 * \author Thomas M. Evans
 * \date   Wed May 07 15:24:31 2014
 * \brief  General_Source test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>

#include "../General_Source.hh"
#include "../Box_Shape.hh"

#include "gtest/utils_gtest.hh"

#include "SourceTestBase.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class GeneralSourceTest : public SourceTestBase
{
    typedef SourceTestBase Base;

  protected:
    typedef profugus::General_Source<profugus::Core> Source;
    typedef Source::SP_Particle                      SP_Particle;
    typedef std::shared_ptr<profugus::Shape>         SP_Shape;

    virtual int get_seed() const
    {
        return 3421;
    }

    virtual void init_group_bounds()
    {
        Vec_Dbl n_bounds;
        n_bounds.push_back(100.);
        n_bounds.push_back(1.);
        n_bounds.push_back(0.001);

        b_group_bounds = std::make_shared<profugus::Group_Bounds>(n_bounds);
    }

    /*
      Test Geometry:

      |-------|-------|
      |       |       |
      |  H20  |  H2O  |
      |       |       |
      |-------|-------|
      |       |       |
      |  H2O  |  UO2  |
      |       |       |
      |-------|-------|

     */
    virtual void init_geometry()
    {
        typedef Geometry_t::SP_Array SP_Core;
        typedef Geometry_t::Array_t  Core_t;
        typedef Core_t::SP_Object    SP_Lattice;
        typedef Core_t::Object_t     Lattice_t;
        typedef Lattice_t::SP_Object SP_Pin_Cell;
        typedef Lattice_t::Object_t  Pin_Cell_t;

        // make cells
        SP_Pin_Cell h2o(std::make_shared<Pin_Cell_t>(3, 1.26, 14.28));

        // make cells
        SP_Pin_Cell uo2(std::make_shared<Pin_Cell_t>(3, 0.54, 3, 1.26, 14.28));

        // make lattice
        SP_Lattice lat(std::make_shared<Lattice_t>(2, 2, 1, 2));
        lat->id(0, 0, 0) = 0;
        lat->id(1, 0, 0) = 1;
        lat->id(0, 1, 0) = 0;
        lat->id(1, 1, 0) = 0;
        lat->assign_object(h2o, 0);
        lat->assign_object(uo2, 1);

        // complete lattice
        lat->complete(0.0, 0.0, 0.0);

        // make the core
        SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));
        core->assign_object(lat, 0);
        core->complete(0.0, 0.0, 0.0);

        // make the b_geometry
        b_geometry = std::make_shared<Geometry_t>(core);
    }

    /*
     - Mat 2 -> UO2
     - Mat 3 -> H20
     */
    virtual void init_physics()
    {
        RCP_XS xs(Teuchos::rcp(new XS_t()));
        xs->set(0, 2);

        XS_t::OneDArray total2(2, 0.5);
        XS_t::TwoDArray scat2(2, 2, 0.1);

        XS_t::OneDArray total3(2, 1.1);
        XS_t::TwoDArray scat3(2, 2, 0.4);
        XS_t::OneDArray bounds(b_group_bounds->group_bounds());

        xs->add(2, XS_t::TOTAL, total2);
        xs->add(2, 0, scat2);
        xs->add(3, XS_t::TOTAL, total3);
        xs->add(3, 0, scat3);
        xs->set_bounds(bounds);

        xs->complete();

        b_physics = std::make_shared<Physics_t>(b_db, xs);
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(GeneralSourceTest, build_and_run)
{
    // Need to run quite a few particles, we're just testing statistics
    // on the results.
    int Np = 200000;
    b_db->set("Np", Np);

    // make a uniform source
    Source source(b_db, b_geometry, b_physics, b_rcon);

    int num_cells  = b_geometry->num_cells();
    int num_groups = b_physics->num_groups();
    std::vector<std::vector<double>> src_vals = {{0.4, 0.2},
                                                 {0.1, 0.2},
                                                 {0.8, 0.4},
                                                 {0.3, 0.3},
                                                 {0.4, 0.6}};

    // build the source
    Teuchos::TwoDArray<double> src(num_cells,num_groups);
    CHECK( src.getNumRows() == 5 );
    CHECK( src.getNumCols() == 2 );

    // Determine expected cell and energy distributions
    std::vector<double> cell_pdf(num_cells,0.0);
    std::vector<std::vector<double>> erg_pdfs(num_cells,
        std::vector<double>(num_groups,0.0));
    for( int cell = 0; cell < num_cells; ++cell )
    {
        auto &erg_pdf = erg_pdfs[cell];
        for( int g = 0; g < num_groups; ++g )
        {
            double val = src_vals[cell][g];
            cell_pdf[cell] += val;
            erg_pdf[g] += val;
            src[cell][g] = val;
        }
        // Norm erg PDF
        double erg_sum = std::accumulate( erg_pdf.begin(), erg_pdf.end(), 0.0 );
        for( auto &val : erg_pdf )
            val /= erg_sum;
    }
    // Normalize cell PDF
    double cell_sum = std::accumulate( cell_pdf.begin(), cell_pdf.end(), 0.0 );
    for( auto &val : cell_pdf )
        val /= cell_sum;

    source.build_source(src);
    EXPECT_TRUE(!source.empty());

    const int nodes = profugus::nodes();
    EXPECT_EQ(Np / nodes, source.num_to_transport());
    EXPECT_EQ(Np, source.total_num_to_transport());
    EXPECT_EQ(0, source.num_run());
    EXPECT_EQ(nodes, source.num_streams());

    std::vector<std::vector<double>> result(num_cells,
        std::vector<double>(num_groups,0.0));
    int ctr = 0;
    while (!source.empty())
    {
        // get a particle
        SP_Particle p = source.get_particle();

        ctr++;

        EXPECT_TRUE(p->alive());
        EXPECT_EQ(1.0, p->wt());

        // Get cell and group info
        int g = p->group();
        int cell = b_geometry->cell(p->geo_state());
        result[cell][g] += 1.0;
    }

    EXPECT_EQ(source.num_run(), ctr);
    EXPECT_EQ(source.num_to_transport(), ctr);
    EXPECT_EQ(source.num_left(), 0);

    // Global reduction of result
    for( int cell = 0; cell < num_cells; ++cell )
    {
        auto &cell_result = result[cell];
        profugus::global_sum(&cell_result[0],num_groups);
    }

    std::vector<double> cell_dist(num_cells,0.0);
    std::vector<std::vector<double>> erg_dists(num_cells,
        std::vector<double>(num_groups,0.0));
    for( int cell = 0; cell < num_cells; ++cell )
    {
        auto &erg_dist = erg_dists[cell];
        for( int g = 0; g < num_groups; ++g )
        {
            double val = result[cell][g];
            cell_dist[cell] += val;
            erg_dist[g] += val;
        }
        // Norm erg PDF
        double erg_sum = std::accumulate( erg_dist.begin(), erg_dist.end(), 0.0 );
        for( auto &val : erg_dist )
            val /= erg_sum;
    }
    // Normalize cell PDF
    cell_sum = std::accumulate( cell_dist.begin(), cell_dist.end(), 0.0 );
    for( auto &val : cell_dist )
        val /= cell_sum;

    double tol = 1e-2;
    EXPECT_VEC_SOFTEQ( cell_pdf, cell_dist, tol );
    for( int cell = 0; cell < num_cells; ++cell )
        EXPECT_VEC_SOFTEQ( erg_pdfs[cell], erg_dists[cell], tol );
}

//---------------------------------------------------------------------------//
//                 end of tstGeneral_Source.cc
//---------------------------------------------------------------------------//
