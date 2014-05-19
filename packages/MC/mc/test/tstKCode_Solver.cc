//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstKCode_Solver.cc
 * \author Thomas M. Evans
 * \date   Mon May 19 13:30:39 2014
 * \brief  KCode_Solver unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../KCode_Solver.hh"

#include "gtest/utils_gtest.hh"

#include <string>
#include <vector>
#include <memory>

#include "harness/DBC.hh"
#include "comm/P_Stream.hh"
#include "xs/XS_Builder.hh"
#include "../Fission_Source.hh"
#include "../Global_RNG.hh"
#include "../Tally.hh"
#include "../Group_Bounds.hh"
#include "../VR_Roulette.hh"

//---------------------------------------------------------------------------//
// Helpers
//---------------------------------------------------------------------------//

class Dummy_Tally : public profugus::Tally
{
    typedef profugus::Tally     Base;
    typedef std::vector<double> Vec_Dbl;

  private:
    int    d_pl_counter;
    double d_finalized_np;

  public:
    // >>> CONSTRUCTION

    // Constructor
    Dummy_Tally(SP_Physics physics)
        : Base(physics)
        , d_pl_counter(1)
        , d_finalized_np(-1.)
    {
        Ensure (d_pl_counter == 1);
    }

    //! Virtual destructor for polymorphism
    virtual ~Dummy_Tally() override final { /* * */ }

    // >>> ACCESSORS

    //! Return the path length counter
    int pl_counter() const { return d_pl_counter; }

    //! Return the value we got called for "finalize"
    int finalized_np() const { return d_finalized_np; }

    // >>> DERIVED INTERFACE

    //! Track particle, using pre-calculated physics information (multipliers)
    inline void accumulate(double step, const Particle_t& p)
    {
        d_pl_counter += 1;
    }

    //! Nothing done at end of particle history
    void end_history() { /* * */ }

    //! Finalize is called at end of program, not cycle, so this is a null-op
    void finalize(double np) { d_finalized_np = np; }

    //! Reset the counter
    void reset() { d_pl_counter = 0; }
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class KCode_SolverTest : public testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::KCode_Solver              Solver_t;
    typedef Solver_t::Physics_t                 Physics_t;
    typedef Solver_t::Geometry_t                Geometry_t;
    typedef profugus::Global_RNG::RNG_Control_t RNG_Control_t;
    typedef Physics_t::Particle_t               Particle_t;
    typedef Physics_t::Bank_t                   Bank_t;
    typedef Solver_t::Tallier_t                 Tallier_t;
    typedef Solver_t::Source_Transporter_t      Transporter_t;
    typedef profugus::VR_Roulette               Var_Reduction_t;

    typedef Physics_t::XS_t   XS_t;
    typedef Physics_t::RCP_XS RCP_XS;

    typedef XS_t::OneDArray OneDArray;
    typedef XS_t::TwoDArray TwoDArray;

    typedef Physics_t::ParameterList_t ParameterList_t;
    typedef Physics_t::RCP_Std_DB      RCP_Std_DB;

    typedef std::shared_ptr<Physics_t>       SP_Physics;
    typedef std::shared_ptr<Geometry_t>      SP_Geometry;
    typedef Physics_t::SP_Particle           SP_Particle;
    typedef std::shared_ptr<RNG_Control_t>   SP_RNG_Control;
    typedef Solver_t::SP_Tallier             SP_Tallier;
    typedef Solver_t::SP_Fission_Source      SP_Fission_Source;
    typedef Solver_t::SP_Source_Transporter  SP_Transporter;
    typedef std::shared_ptr<Var_Reduction_t> SP_Var_Reduction;

  protected:
    void SetUp()
    {
        seed = 32442;
        rcon = std::make_shared<RNG_Control_t>(seed);

        // set number of nodes
        node  = profugus::node();
        nodes = profugus::nodes();

        // set our other attributes
        this->init_db();
        this->init_geometry();
        this->init_physics();
    }

    //! Set the geometry
    void init_geometry()
    {
        typedef Geometry_t::SP_Array SP_Core;
        typedef Geometry_t::Array_t  Core_t;
        typedef Core_t::SP_Object    SP_Lattice;
        typedef Core_t::Object_t     Lattice_t;
        typedef Lattice_t::SP_Object SP_Pin_Cell;
        typedef Lattice_t::Object_t  Pin_Cell_t;

        // UO2 pin cell
        SP_Pin_Cell pin(std::make_shared<Pin_Cell_t>(1, 0.54, 8, 1.26, 1.0));

        // make lattice
        SP_Lattice lat(std::make_shared<Lattice_t>(1, 1, 1, 1));

        // assign pins
        lat->assign_object(pin, 0); // fuel pins

        // arrange pin-cell in lattice
        lat->id(0, 0, 0) = 0; // fuel pin

        // complete lattice
        lat->complete(0.0, 0.0, 0.0);

        // make core
        SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));

        // assign lattice
        core->assign_object(lat, 0); // lattice of single fuel pin

        // arrange lattice in core
        core->id(0, 0, 0) = 0;

        // add reflecting boundary conditions
        Core_t::Vec_Int reflect(6, 1);
        core->set_reflecting(reflect);

        // complete core
        core->complete(0.0, 0.0, 0.0);

        // make the geometry
        geometry = std::make_shared<Geometry_t>(core);
    }

    //! Set the physics
    void init_physics()
    {
        typedef profugus::XS_Builder::Matid_Map Matid_Map;
        typedef Matid_Map::value_type           value_type;
        typedef std::string                     string;

        // make an xs_builder
        profugus::XS_Builder builder;
        builder.open_and_broadcast("xs_c5g7.xml");

        Matid_Map map;
        map.insert(value_type(1, string("UO2")));
        map.insert(value_type(2, string("MOX_43")));
        map.insert(value_type(3, string("MOX_70")));
        map.insert(value_type(4, string("MOX_87")));
        map.insert(value_type(5, string("fission_chamber")));
        map.insert(value_type(6, string("guide_tube")));
        map.insert(value_type(7, string("control_rod")));
        map.insert(value_type(8, string("mod")));
        map.complete();

        builder.build(map);

        RCP_XS xs = builder.get_xs();

        physics = std::make_shared<Physics_t>(db, xs);
        physics->set_geometry(geometry);
    }

    //! Set the database options.
    void init_db()
    {
        db = Teuchos::rcp(new ParameterList_t("test"));

        db->set("Np", 1000);
        OneDArray fs(6, 0.09);
        fs[1] = 1.17; fs[3] = 1.17; fs[4] = 0.0; fs[5] = 1.00;
        db->set("init_fission_src", fs);

        db->set("mc_diag_frac", 1.1);
        db->set("keff_init", 1.32);
        db->set("num_cycles", 20);
        db->set("num_inactive_cycles", 10);
    }

  protected:
    int node;
    int nodes;
    int seed;

    SP_RNG_Control  rcon;
    RCP_Std_DB      db;
    SP_Geometry     geometry;
    SP_Physics      physics;

    SP_Transporter   transporter;
    SP_Var_Reduction var_reduction;
    SP_Tallier       tallier;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(KCode_SolverTest, pin_cell)
{
    // make the source transporter
    transporter = std::make_shared<Transporter_t>(db, geometry, physics);

    // make the variance reduction
    var_reduction = std::make_shared<Var_Reduction_t>(db);

    // make the tallier
    tallier = std::make_shared<Tallier_t>();
    tallier->set(geometry, physics);;

    // add objects to the source transporter
    transporter->set(tallier);
    transporter->set(var_reduction);

    // make Kcode-solver
    Solver_t solver(db);

    // make fission source
    SP_Fission_Source fsrc(std::make_shared<profugus::Fission_Source>(
                               db, geometry, physics, rcon));

    // set it
    solver.set(transporter, fsrc);

    // solve
    solver.solve();

    // get keff output from tally
    const auto &keff_tally = *solver.keff_tally();

    EXPECT_EQ(10, keff_tally.cycle_count());
    auto keff = keff_tally.latest();
    profugus::global_barrier();

    profugus::pcout << "Final Keff = " << profugus::setw(12)
                    << profugus::fixed << keff
                    << profugus::endl;
    profugus::pcout << "Note: this is an uncoverged result."
                    << profugus::endl << profugus::endl;

    profugus::global_barrier();

    EXPECT_SOFTEQ(1.3, keff, 0.10);
    if (node == 0)
        EXPECT_SOFTEQ(1.3, keff_tally.mean(), 0.05);
}

//---------------------------------------------------------------------------//
//                 end of tstKCode_Solver.cc
//---------------------------------------------------------------------------//
