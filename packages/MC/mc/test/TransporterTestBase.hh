//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/TransporterTestBase.hh
 * \author Seth R Johnson and Tom Evans
 * \date   Monday May 12 14:20:38 2014
 * \brief  TransporterTestBase class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_test_TransporterTestBase_hh
#define mc_test_TransporterTestBase_hh

#include <Utils/config.h>
#include "gtest/gtest.h"

#include <vector>
#include <memory>

#include "../Shape.hh"
#include "../Global_RNG.hh"
#include "../Group_Bounds.hh"
#include "../Physics.hh"
#include "../Tallier.hh"
#include "../Source.hh"
#include "../Variance_Reduction.hh"
#include "geometry/RTK_Geometry.hh"

//===========================================================================//
/*!
 * \class TransporterTestBase
 * \brief Base class for Shift transport testing.
 */
//===========================================================================//

class TransporterTestBase : public testing::Test
{
    typedef testing::Test Base;

  protected:
    typedef profugus::Physics<profugus::Core>   Physics_t;
    typedef Physics_t::Geometry_t               Geometry_t;
    typedef profugus::Global_RNG::RNG_Control_t RNG_Control_t;
    typedef Physics_t::Particle_t               Particle_t;
    typedef Physics_t::SP_Particle              SP_Particle;
    typedef Physics_t::Bank_t                   Bank_t;

    typedef Physics_t::XS_t   XS_t;
    typedef Physics_t::RCP_XS RCP_XS;

    typedef Physics_t::ParameterList_t ParameterList_t;
    typedef Physics_t::RCP_Std_DB      RCP_Std_DB;

    typedef Physics_t::Space_Vector Space_Vector;
    typedef Physics_t::Space_Vector Vector;

    typedef std::vector<double> Vec_Dbl;

    typedef profugus::Tallier            Tallier_t;
    typedef profugus::Variance_Reduction Var_Reduction_t;
    typedef profugus::Source<Geometry_t> Source_t;

    typedef std::shared_ptr<Geometry_t>             SP_Geometry;
    typedef std::shared_ptr<Physics_t>              SP_Physics;
    typedef std::shared_ptr<RNG_Control_t>          SP_RNG_Control;
    typedef std::shared_ptr<profugus::Shape>        SP_Shape;
    typedef std::shared_ptr<Tallier_t>              SP_Tallier;
    typedef std::shared_ptr<Var_Reduction_t>        SP_Var_Reduction;
    typedef std::shared_ptr<profugus::Group_Bounds> SP_Group_Bounds;

  protected:
    // Initialization that are performed for each test
    virtual void SetUp();

    // 1 group
    virtual void init_group_bounds();

    // Create group information
    virtual void init_db();

    /* 4 material definitions/1 group
     *
     - Mat 0 -> Undefined
     - Mat 1 -> Fuel
     - Mat 2 -> Guide tube
     - Mat 3 -> Moderater
     */
    virtual void init_geometry();

    // Set physics with total and self-scattering
    virtual void init_physics();

    // Create russian roulette VR by default
    virtual void init_vr();

    // Set tallies
    virtual void init_tallies()
    {
        /* default is to have no tallies */
    }

    // Number of groups
    int num_groups() const
    {
        return group_bounds->num_groups();
    }

  protected:
    int node;
    int nodes;

    SP_RNG_Control   rcon;
    RCP_Std_DB       db;
    SP_Geometry      geometry;
    SP_Physics       physics;
    SP_Group_Bounds  group_bounds;
    SP_Var_Reduction var_red;
    SP_Tallier       tallier;
};

#endif // mc_test_TransporterTestBase_hh

//---------------------------------------------------------------------------//
//              end of mc_transport/TransporterTestBase.hh
//---------------------------------------------------------------------------//
