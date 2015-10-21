//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstFixed_Source_Solver.cc
 * \author Thomas M. Evans
 * \date   Tue May 13 17:11:19 2014
 * \brief  Fixed_Source_Solver test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Fixed_Source_Solver.hh"

#include "gtest/utils_gtest.hh"

#include <memory>

#include "../Box_Shape.hh"
#include "../Uniform_Source.hh"
#include "../VR_Roulette.hh"

#include "TransporterTestBase.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class FixedSourceSolverTest : public TransporterTestBase
{
    typedef TransporterTestBase Base;

  public:
    typedef profugus::Source_Transporter  Transporter_t;
    typedef profugus::Fixed_Source_Solver Solver_t;
    typedef Transporter_t::Source_t       Source_t;
    typedef Solver_t::SP_Source           SP_Source;

    typedef std::shared_ptr<Transporter_t> SP_Transporter;

  protected:
    void SetUp()
    {
        Base::SetUp();

        db->set("Np", 10000);
        db->set("mc_diag_frac", 0.2);

        this->init_source();

        transporter = std::make_shared<Transporter_t>(db, geometry, physics);
        transporter->set(tallier);
        transporter->set(var_red);
    }

    void init_geometry()
    {
        typedef Geometry_t::SP_Array SP_Core;
        typedef Geometry_t::Array_t  Core_t;
        typedef Core_t::SP_Object    SP_Lattice;
        typedef Core_t::Object_t     Lattice_t;
        typedef Lattice_t::SP_Object SP_Pin_Cell;
        typedef Lattice_t::Object_t  Pin_Cell_t;

        // box
        SP_Pin_Cell box(std::make_shared<Pin_Cell_t>(0, 5.0, 5.0));

        // make lattice
        SP_Lattice lat(std::make_shared<Lattice_t>(1, 1, 1, 1));

        // assign pins
        lat->assign_object(box, 0); // uniform lattice

        // complete lattice
        lat->complete(0.0, 0.0, 0.0);

        // make core
        SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));
        core->assign_object(lat, 0);

        // add reflecting boundary conditions
        Core_t::Vec_Int reflect(6, 1);
        core->set_reflecting(reflect);

        core->complete(0.0, 0.0, 0.0);

        geometry = std::make_shared<Geometry_t>(core);
    }

    void init_group_bounds()
    {
        Vec_Dbl n_bounds;
        n_bounds.push_back(2.e7);
        n_bounds.push_back(1.e-5);

        group_bounds = std::make_shared<profugus::Group_Bounds>(n_bounds);
    }

    void init_physics()
    {
        RCP_XS xs(Teuchos::rcp(new XS_t()));
        xs->set(0, 1);

        XS_t::OneDArray tot1(1, 5.0);
        XS_t::TwoDArray sct1(1, 1, 1.5);
        XS_t::OneDArray bounds(group_bounds->group_bounds());

        xs->set_bounds(bounds);
        xs->add(0, XS_t::TOTAL, tot1);
        xs->add(0, 0, sct1);

        xs->complete();

        physics = std::make_shared<Physics_t>(db, xs);
    }

    void init_source()
    {
        auto box = std::make_shared<profugus::Box_Shape>(
            0.0, 5.0, 0.0, 5.0, 0.0, 5.0);

        auto usource = std::make_shared<profugus::Uniform_Source>(
            db, geometry, physics, box );

        source = usource;
    }

    void init_vr()
    {
        var_red = std::make_shared<profugus::VR_Roulette>(db);
    }

  protected:
    SP_Source      source;
    SP_Transporter transporter;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(FixedSourceSolverTest, inf_med)
{
#if 0 // NO TALLIES YET
       ASSERT_TRUE(cell_tally);
       ASSERT_EQ(1, tallier->num_pathlength_tallies());
#endif

    // we may not divide evenly depending on number of nodes
    ASSERT_NEAR(10000, source->total_num_to_transport(), 10);

    // make the fixed solver
    Solver_t solver;

    solver.set(transporter, source);
    solver.solve();

#if 0 // NO TALLIES YET
    // ref = Q / sig_a = 1.0/cc / 3.5/cc
    double ref = 1.0 / 3.5;
    double src = 1.0 * 125.0;

    // we don't normalize by volume yet, but the MC result is per source
    // neutron
    double ans = cell_tally->result()->mean(0)[0] / 125.0;
    double flux = ans * src;

    double var = cell_tally->result()->variance(0)[0] * src / 125.;
    double re = std::sqrt(var) / flux;

    EXPECT_SOFTEQ(ref, flux, re * 2); // 95% confidence
    EXPECT_SOFTEQ(ref, flux, 0.016);  // assuming ~10k particles
#endif
}

//---------------------------------------------------------------------------//
//                 end of tstFixed_Source_Solver.cc
//---------------------------------------------------------------------------//
