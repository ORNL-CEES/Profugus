//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/SourceTestBase.hh
 * \author Seth R Johnson and Tom Evans
 * \date   Tuesday May 6 12:0:29 2014
 * \brief  Source test infrastructure base class.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_test_SourceTestBase_hh
#define mc_test_SourceTestBase_hh

#include <Utils/config.h>
#include "gtest/gtest.h"

#include <vector>
#include <memory>

#include "harness/DBC.hh"
#include "rng/RNG.hh"
#include "../Group_Bounds.hh"
#include "../Physics.hh"

//===========================================================================//
/*!
 * \class SourceTestBase
 * \brief Base class for google test fixture
 *
 * Any of the construction methods can be overridden by a subclass.
 */
//===========================================================================//

class SourceTestBase : public testing::Test
{
  protected:
    typedef profugus::Physics                   Physics_t;
    typedef Physics_t::Geometry_t               Geometry_t;
    typedef profugus::RNG                       RNG_t;

    typedef Physics_t::XS_t   XS_t;
    typedef Physics_t::RCP_XS RCP_XS;

    typedef Physics_t::ParameterList_t ParameterList_t;
    typedef Physics_t::RCP_Std_DB      RCP_Std_DB;

    typedef std::shared_ptr<Geometry_t>    SP_Geometry;
    typedef std::shared_ptr<Physics_t>     SP_Physics;
    typedef std::shared_ptr<RNG_t> SP_RNG;

    typedef std::shared_ptr<profugus::Group_Bounds> SP_Group_Bounds;

    typedef Physics_t::Space_Vector Space_Vector;

    typedef std::vector<double> Vec_Dbl;

  protected:
    //! Initialization that are performed for each test
    virtual void SetUp();

    // Get random number seed
    virtual int get_seed() const;

    // Construct MG physics energy boundaries
    virtual void init_group_bounds();

    // Create group information
    virtual void init_db();

    // Fill with four vertical cells of water from -2,-2,-4 to 2,2,4
    virtual void init_geometry();

    // Set physics with total and self-scattering
    virtual void init_physics();

    // Number of groups
    int num_groups() const
    {
        return b_group_bounds->num_groups();
    }

  protected:
    int             node;
    int             nodes;
    SP_RNG          b_rcon;
    RCP_Std_DB      b_db;
    SP_Geometry     b_geometry;
    SP_Physics      b_physics;
    SP_Group_Bounds b_group_bounds;
};

#endif // mc_test_SourceTestBase_hh

//---------------------------------------------------------------------------//
//              end of mc_sources/test/SourceTestBase.hh
//---------------------------------------------------------------------------//
