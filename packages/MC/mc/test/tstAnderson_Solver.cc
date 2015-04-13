//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/test/tstAnderson_Solver.cc
 * \author Thomas M. Evans
 * \date   Mon Apr 13 14:00:16 2015
 * \brief  Anderson_Solver class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Anderson_Solver.hh"

#include "gtest/utils_gtest.hh"

#include <string>
#include <vector>
#include <memory>

#include "xs/XS_Builder.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "spn/MatrixTraits.hh"
#include "geometry/Cartesian_Mesh.hh"
#include "../Fission_Source.hh"
#include "../Global_RNG.hh"
#include "../Tally.hh"
#include "../Group_Bounds.hh"
#include "../VR_Roulette.hh"
#include "../Tallier.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class Anderson_SolverTest : public testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::EpetraTypes               LinAlg_t;
    typedef profugus::Anderson_Solver<LinAlg_t> Solver_t;
    typedef Solver_t::Operator                  Operator;
    typedef Operator::Source_Transporter_t      Transporter_t;
    typedef profugus::Tallier                   Tallier_t;
    typedef Transporter_t::Physics_t            Physics_t;
    typedef Transporter_t::Geometry_t           Geometry_t;
    typedef profugus::Global_RNG::RNG_Control_t RNG_Control_t;
    typedef Physics_t::Particle_t               Particle_t;
    typedef Physics_t::Bank_t                   Bank_t;
    typedef profugus::VR_Roulette               Var_Reduction_t;
    typedef profugus::Cartesian_Mesh            Mesh_t;

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

    typedef Solver_t::RCP_MAP RCP_MAP;

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

        // make the source transporter
        transporter = std::make_shared<Transporter_t>(db, geometry, physics);

        // make the variance reduction
        var_reduction = std::make_shared<Var_Reduction_t>(db);

        // make the tallier
        tallier = std::make_shared<Tallier_t>();
        tallier->set(geometry, physics);

        // add objects to the source transporter
        transporter->set(tallier);
        transporter->set(var_reduction);
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
        SP_Lattice lat(std::make_shared<Lattice_t>(3, 3, 1, 1));

        // assign pins
        lat->assign_object(pin, 0); // fuel pins

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

        db->set("mc_diag_frac", 1.1);
        db->set("keff_init", 1.32);
        db->set("num_cycles", 20);

        auto &adb = db->sublist("anderson_db");
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

TEST_F(Anderson_SolverTest, api)
{
    // Make anderson mesh
    std::vector<double> r = {0.0, 1.26, 2.52, 3.78};
    std::vector<double> z = {0.0, 1.0};
    auto &adb = db->sublist("anderson_db");
    adb.set("x_bounds", OneDArray(r));
    adb.set("y_bounds", OneDArray(r));
    adb.set("z_bounds", OneDArray(z));
    
    // make fission source
    auto fs = std::make_shared<profugus::Fission_Source>(
        db, geometry, physics, rcon);

    // Make the operator
    Solver_t anderson(db);

    // Set it
    anderson.set(transporter, fs);
}

//---------------------------------------------------------------------------//
// end of MC/mc/test/tstAnderson_Solver.cc
//---------------------------------------------------------------------------//
