//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/SourceTestBase.cc
 * \author Seth R Johnson
 * \date   Tuesday May 6 12:0:22 2014
 * \brief  Source test infrastructure class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "SourceTestBase.hh"

#include "Teuchos_RCP.hpp"

#include "comm/global.hh"

//---------------------------------------------------------------------------//

//! Initialization that are performed for each test
void SourceTestBase::SetUp()
{
    b_rcon = std::make_shared<RNG_Control_t>(this->get_seed());

    // set number of nodes
    node  = profugus::node();
    nodes = profugus::nodes();

    // set our other attributes
    this->init_group_bounds();
    this->init_db();
    this->init_geometry();
    this->init_physics();

    CHECK(b_physics);
    CHECK(b_geometry);
    b_physics->set_geometry(b_geometry);

    ENSURE(b_rcon);
    ENSURE(!b_db.is_null());
    ENSURE(b_geometry);
    ENSURE(b_physics);
    ENSURE(b_group_bounds);
}

//---------------------------------------------------------------------------//

//! Get random number seed
int SourceTestBase::get_seed() const
{
    return 48151623;
}

//---------------------------------------------------------------------------//

//! Construct energy boundaries
void SourceTestBase::init_group_bounds()
{
    Vec_Dbl n_bounds;
    n_bounds.push_back(2.e7);
    n_bounds.push_back(1.e-3);

    b_group_bounds = std::make_shared<profugus::Group_Bounds>(n_bounds);

    ENSURE(b_group_bounds->num_groups() == 1);
}

//---------------------------------------------------------------------------//

//! Create group information
void SourceTestBase::init_db()
{
    REQUIRE(b_db.is_null());
    REQUIRE(b_group_bounds);

    b_db = Teuchos::rcp(new ParameterList_t("test"));
}

//---------------------------------------------------------------------------//

//! Fill with four vertical cells of water from -2,-2,-4 to 2,2,4
void SourceTestBase::init_geometry()
{
    REQUIRE(!b_db.is_null());
    REQUIRE(!b_geometry);

    typedef Geometry_t::SP_Array SP_Core;
    typedef Geometry_t::Array_t  Core_t;
    typedef Core_t::SP_Object    SP_Lattice;
    typedef Core_t::Object_t     Lattice_t;
    typedef Lattice_t::SP_Object SP_Pin_Cell;
    typedef Lattice_t::Object_t  Pin_Cell_t;

    // make cells
    // (material pitch height)
    SP_Pin_Cell h2o(std::make_shared<Pin_Cell_t>(3, 2., 8.));

    // make lattice
    // (nx ny nz num_objects)
    SP_Lattice lat(std::make_shared<Lattice_t>(2, 2, 1, 1));
    lat->assign_object(h2o, 0);

    // complete lattice
    lat->complete(0.0, 0.0, 0.0);

    // make core
    SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));
    core->assign_object(lat, 0);
    core->complete(-2., -2., -4.); //lower corner is -2,-2,-4

    // make the b_geometry
    b_geometry = std::make_shared<Geometry_t>(core);
}

//---------------------------------------------------------------------------//

//! Set materials with scattering
void SourceTestBase::init_physics()
{
    REQUIRE(!b_db.is_null());
    REQUIRE(b_group_bounds);
    const int ng = num_groups();

    RCP_XS xs(Teuchos::rcp(new XS_t()));
    xs->set(0, ng);

    const double tot_scattering  = 0.9;
    const double self_scattering = 0.8;

    XS_t::OneDArray total(ng, 0.0);
    XS_t::TwoDArray scat(ng, ng, 0.0);

    for (int g = 0; g < ng; ++g)
    {
        total[g]   = g + 0.1;
        scat(g, g) = total[g] * self_scattering;

        if (g > 0)
        {
            scat(g, g-1) = total[g] * (tot_scattering - self_scattering);
        }
    }

    XS_t::OneDArray bounds(b_group_bounds->group_bounds());

    xs->set_bounds(bounds);
    xs->add(3, XS_t::TOTAL, total);
    xs->add(3, 0, scat);

    xs->complete();

    b_physics = std::make_shared<Physics_t>(b_db, xs);
}

//---------------------------------------------------------------------------//
//              end of mc_sources/test/SourceTestBase.cc
//---------------------------------------------------------------------------//
