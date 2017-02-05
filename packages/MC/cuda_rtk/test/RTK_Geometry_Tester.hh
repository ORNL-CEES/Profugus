//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_rtk/test/RTK_Geometry_Tester.hh
 * \author Tom Evans
 * \date   Fri Feb 03 09:50:55 2017
 * \brief  RTK_Geometry_Tester class declaration.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_test_RTK_Geometry_Tester_hh
#define MC_cuda_rtk_test_RTK_Geometry_Tester_hh

#include <vector>
#include <random>
#include <limits>
#include <memory>

#include "Utils/gtest/Gtest_Functions.hh"
#include "MC/geometry/RTK_Geometry.hh"

//---------------------------------------------------------------------------//
// TESTERS
//---------------------------------------------------------------------------//

class Base : public ::testing::Test
{
  protected:
    // Types
    using Core_Geometry = profugus::Core;
    using SP_Geometry   = std::shared_ptr<Core_Geometry>;
    using Core_t        = Core_Geometry::Array_t;
    using Lattice_t     = Core_t::Object_t;
    using Pin_Cell_t    = Lattice_t::Object_t;
    using SP_Core       = Core_Geometry::SP_Array;
    using SP_Lattice    = Core_t::SP_Object;
    using SP_Pin_Cell   = Lattice_t::SP_Object;

  protected:

    virtual ~Base() = default;

  protected:

    SP_Geometry geometry;
};

//---------------------------------------------------------------------------//

class Core : public Base
{
  protected:

    void SetUp()
    {
        // 2 fuel pin types
        SP_Pin_Cell pin1(std::make_shared<Pin_Cell_t>(1, 0.54, 3, 1.26, 14.28));
        SP_Pin_Cell pin2(std::make_shared<Pin_Cell_t>(2, 0.54, 3, 1.26, 14.28));

        // water pin
        SP_Pin_Cell box(std::make_shared<Pin_Cell_t>(3, 2.52, 14.28));

        // 3 lattices (fuel 1, 2, and water)
        SP_Lattice lat1(std::make_shared<Lattice_t>(2, 2, 1, 1));
        SP_Lattice lat2(std::make_shared<Lattice_t>(2, 2, 1, 1));
        SP_Lattice lat3(std::make_shared<Lattice_t>(1, 1, 1, 1));

        // lattice assignments
        lat1->assign_object(pin1, 0);
        lat2->assign_object(pin2, 0);
        lat3->assign_object(box, 0);

        // complete the lattices
        lat1->complete(0.0, 0.0, 0.0);
        lat2->complete(0.0, 0.0, 0.0);
        lat3->complete(0.0, 0.0, 0.0);

        // make core (3x3x2 with 4 objects, object 0 unassigned)
        SP_Core core(std::make_shared<Core_t>(3, 3, 2, 4));
        EXPECT_EQ(1, core->level());

        // assign lattices
        core->assign_object(lat1, 1);
        core->assign_object(lat2, 2);
        core->assign_object(lat3, 3);

        // assign ids
        core->id(0, 0, 0) = 1; // lattice 1
        core->id(1, 0, 0) = 2; // lattice 2
        core->id(2, 0, 0) = 3; // lattice 3 (reflector)
        core->id(0, 1, 0) = 2; // lattice 2
        core->id(1, 1, 0) = 1; // lattice 1
        core->id(2, 1, 0) = 3; // lattice 3 (reflector)
        core->id(0, 2, 0) = 3; // lattice 3 (reflector)
        core->id(1, 2, 0) = 3; // lattice 3 (reflector)
        core->id(2, 2, 0) = 3; // lattice 3 (reflector)
        core->id(0, 0, 1) = 3; // lattice 3 (reflector)
        core->id(1, 0, 1) = 3; // lattice 3 (reflector)
        core->id(2, 0, 1) = 3; // lattice 3 (reflector)
        core->id(0, 1, 1) = 3; // lattice 3 (reflector)
        core->id(1, 1, 1) = 3; // lattice 3 (reflector)
        core->id(2, 1, 1) = 3; // lattice 3 (reflector)
        core->id(0, 2, 1) = 3; // lattice 3 (reflector)
        core->id(1, 2, 1) = 3; // lattice 3 (reflector)
        core->id(2, 2, 1) = 3; // lattice 3 (reflector)

        // complete the core
        core->complete(0.0, 0.0, 0.0);

        // finish the geometry
        geometry = std::make_shared<Core_Geometry>(core);

        // Set number of threads to run on
        num_threads = 256;
        num_blocks  = 8;

        // Make random number seeds
        seeds.resize(num_threads * num_blocks);
        {
            std::mt19937 rng(42342);
            std::uniform_int_distribution<unsigned long> dist(
                std::numeric_limits<unsigned long>::min(),
                std::numeric_limits<unsigned long>::max());

            for (auto &seed : seeds)
            {
                seed = dist(rng);
            }
        }
    }

    void heuristic();

  protected:

    std::vector<unsigned long> seeds;
    unsigned int num_threads;
    unsigned int num_blocks;
};

#endif // MC_cuda_rtk_test_RTK_Geometry_Tester_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/test/RTK_Geometry_Tester.hh
//---------------------------------------------------------------------------//
