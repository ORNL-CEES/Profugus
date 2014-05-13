 //----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/TransporterTestBase.cc
 * \author Seth R Johnson and Thomas M. Evans
 * \date   Fri Mar 29 17:45:16 2013
 * \brief  TransporterTestBase member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "TransporterTestBase.hh"

#include "Teuchos_RCP.hpp"

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"

#include "../VR_Roulette.hh"

//---------------------------------------------------------------------------//
//! Initialization that are performed for each test

void TransporterTestBase::SetUp()
{
    rcon = std::make_shared<profugus::RNG_Control>(349832);

    // set number of nodes
    node  = profugus::node();
    nodes = profugus::nodes();

    this->init_group_bounds();
    Check(group_bounds);

    this->init_db();
    Check(!db.is_null());

    this->init_geometry();
    Check(geometry);

    this->init_physics();
    Check(physics);
    physics->set_geometry(geometry);

    // Create VR
    this->init_vr();
    Check(var_red);
    var_red->set(geometry);
    var_red->set(physics);

    // Create tallier
    tallier = std::make_shared<Tallier_t>();
    tallier->set(geometry, physics);
    this->init_tallies();

    Ensure(rcon);
    Ensure(!db.is_null());
    Ensure(geometry);
    Ensure(physics);
    Ensure(group_bounds);
    Ensure(var_red);
    Ensure(tallier);
}

//---------------------------------------------------------------------------//

void TransporterTestBase::init_group_bounds()
{
    Vec_Dbl n_bounds;
    n_bounds.push_back(100.);
    n_bounds.push_back(0.001);

    group_bounds = std::make_shared<profugus::Group_Bounds>(n_bounds);

    Ensure(group_bounds->num_groups() == 1);
}

//---------------------------------------------------------------------------//

//! Create group information
void TransporterTestBase::init_db()
{
    Require(db.is_null());
    Require(group_bounds);

    db = Teuchos::rcp(new ParameterList_t("test"));
}

//---------------------------------------------------------------------------//

void TransporterTestBase::init_geometry()
{
    typedef Geometry_t::SP_Array SP_Core;
    typedef Geometry_t::Array_t  Core_t;
    typedef Core_t::SP_Object    SP_Lattice;
    typedef Core_t::Object_t     Lattice_t;
    typedef Lattice_t::SP_Object SP_Pin_Cell;
    typedef Lattice_t::Object_t  Pin_Cell_t;

    // make pin cells
    SP_Pin_Cell p1(std::make_shared<Pin_Cell_t>(1, 0.54, 3, 1.26, 14.28));
    SP_Pin_Cell p2(std::make_shared<Pin_Cell_t>(2, 0.54, 3, 1.26, 14.28));

    // make lattice
    SP_Lattice lat(std::make_shared<Lattice_t>(3, 3, 1, 2));

    // assign pins
    lat->assign_object(p1, 0); // fuel pins
    lat->assign_object(p2, 1); // guide tube

    // arrange pin-cells in lattice
    lat->id(0, 0, 0) = 0; // fuel pin
    lat->id(1, 0, 0) = 0; // fuel pin
    lat->id(2, 0, 0) = 0; // fuel pin
    lat->id(0, 1, 0) = 0; // fuel pin
    lat->id(1, 1, 0) = 1; // guide tube
    lat->id(2, 1, 0) = 0; // fuel pin
    lat->id(0, 2, 0) = 0; // fuel pin
    lat->id(1, 2, 0) = 0; // fuel pin
    lat->id(2, 2, 0) = 0; // fuel pin

    // complete lattice
    lat->complete(0.0, 0.0, 0.0);

    Check(profugus::soft_equiv(lat->pitch(def::X), 3.78));
    Check(profugus::soft_equiv(lat->pitch(def::Y), 3.78));
    Check(profugus::soft_equiv(lat->height(), 14.28));

    // make core
    SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));
    core->assign_object(lat, 0);
    core->complete(0.0, 0.0, 0.0);

    geometry = std::make_shared<Geometry_t>(core);
}

//---------------------------------------------------------------------------//

void TransporterTestBase::init_physics()
{
    Require (!db.is_null());
    Require (group_bounds);
    const int ng = num_groups();

    RCP_XS xs(Teuchos::rcp(new XS_t()));
    xs->set(0, ng);

    XS_t::OneDArray tot1(ng, 10.0);
    XS_t::OneDArray tot2(ng, 1.5);
    XS_t::OneDArray tot3(ng, 1.1);

    XS_t::TwoDArray sct1(ng, ng, 4.1);
    XS_t::TwoDArray sct2(ng, ng, 1.4);
    XS_t::TwoDArray sct3(ng, ng, 0.9);

    XS_t::OneDArray bounds(group_bounds->group_bounds());

    xs->set_bounds(bounds);
    xs->add(1, XS_t::TOTAL, tot1);
    xs->add(2, XS_t::TOTAL, tot2);
    xs->add(3, XS_t::TOTAL, tot3);
    xs->add(1, 0, sct1);
    xs->add(2, 0, sct2);
    xs->add(3, 0, sct3);

    XS_t::OneDArray fis(ng, 4.0);
    XS_t::OneDArray chi(ng, 1.0);
    XS_t::OneDArray nuf(ng, 2.4*4.0);

    xs->add(1, XS_t::SIG_F, fis);
    xs->add(1, XS_t::NU_SIG_F, nuf);
    xs->add(1, XS_t::CHI, chi);

    xs->complete();

    physics = std::make_shared<Physics_t>(db, xs);
}

//---------------------------------------------------------------------------//

void TransporterTestBase::init_vr()
{
    db->set("weight_cutoff", 0.01);

    var_red = std::make_shared<profugus::VR_Roulette>(db);
}

//---------------------------------------------------------------------------//
//                 end of TransporterTestBase.cc
//---------------------------------------------------------------------------//
