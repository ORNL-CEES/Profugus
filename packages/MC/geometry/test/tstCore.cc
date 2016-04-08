//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/test/tstRTK_Core.cc
 * \author Thomas M. Evans
 * \date   Mon Jan 24 09:48:47 2011
 * \brief  RTK_Core unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <memory>

#include "utils/Definitions.hh"
#include "utils/Constants.hh"
#include "utils/Vector_Functions.hh"
#include "rng/RNG_Control.hh"
#include "../Definitions.hh"
#include "../RTK_Geometry.hh"

using namespace std;

using def::X;
using def::Y;
using def::Z;

typedef profugus::Core          Core_Geometry;
typedef Core_Geometry           Geometry;
typedef Core_Geometry::Array_t  Core_t;
typedef Core_t::Object_t        Lattice_t;
typedef Lattice_t::Object_t     Pin_Cell_t;
typedef Core_Geometry::SP_Array SP_Core;
typedef Core_t::SP_Object       SP_Lattice;
typedef Lattice_t::SP_Object    SP_Pin_Cell;

typedef Geometry::Space_Vector Vector;
typedef Geometry::Geo_State_t  State;

using profugus::geometry::OUTSIDE;
using profugus::geometry::INSIDE;
using profugus::geometry::REFLECT;

int seed = 4305834;

bool do_output = false;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
/*
 Core Model

 Fuel region         = 1
 Fuel region         = 2
 Moderator/Reflector = 3

 =============================
 ||       ||       ||       ||
 ||       ||       ||       ||
 ||   3   ||   3   ||   3   ||
 ||       ||       ||       ||
 ||       ||       ||       ||
 =============================
 || 2 | 2 || 1 | 1 ||       ||
 ||   |   ||   |   ||       ||
 ||-------||-------||   3   ||
 || 2 | 2 || 1 | 1 ||       ||
 ||   |   ||   |   ||       ||
 =============================
 || 1 | 1 || 2 | 2 ||       ||
 ||   |   ||   |   ||       ||
 ||-------||-------||   3   ||
 || 1 | 1 || 2 | 2 ||       ||
 ||   |   ||   |   ||       ||
 =============================

 */

TEST(Core, Heuristic)
{
    // 2 fuel pin types
    SP_Pin_Cell pin1(make_shared<Pin_Cell_t>(1, 0.54, 3, 1.26, 14.28));
    SP_Pin_Cell pin2(make_shared<Pin_Cell_t>(2, 0.54, 3, 1.26, 14.28));

    // water pin
    SP_Pin_Cell box(make_shared<Pin_Cell_t>(3, 2.52, 14.28));

    // 3 lattices (fuel 1, 2, and water)
    SP_Lattice lat1(make_shared<Lattice_t>(2, 2, 1, 1));
    SP_Lattice lat2(make_shared<Lattice_t>(2, 2, 1, 1));
    SP_Lattice lat3(make_shared<Lattice_t>(1, 1, 1, 1));

    // lattice assignments
    lat1->assign_object(pin1, 0);
    lat2->assign_object(pin2, 0);
    lat3->assign_object(box, 0);

    // complete the lattices
    lat1->complete(0.0, 0.0, 0.0);
    lat2->complete(0.0, 0.0, 0.0);
    lat3->complete(0.0, 0.0, 0.0);

    // make core (3x3x2 with 4 objects, object 0 unassigned)
    SP_Core core(make_shared<Core_t>(3, 3, 2, 4));
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

    // check lattice
    EXPECT_TRUE(soft_equiv(core->pitch(X), 7.56));
    EXPECT_TRUE(soft_equiv(core->pitch(Y), 7.56));
    EXPECT_TRUE(soft_equiv(core->height(), 28.56));

    // output the core array
    core->output(cout);
    cout << endl;

    // build the RTK Core
    Core_Geometry rtk_core(core);
    Geometry &geometry = rtk_core;

    // make a random number generator
    profugus::RNG_Control control(seed);
    auto rng = control.rng();

    // plot collision sites
    ofstream csites("csites.dat");

    // geometry variables
    double costheta, sintheta, phi;
    Vector r, omega;
    State  state;
    double d;
    int    Np = 10000;

    int face_bin[6] = {0};

    // sample Np tracks
    for (int n = 0; n < Np; ++n)
    {
        // sample x,y,z randomly
        r[0] = rng.ran() * 7.56;
        r[1] = rng.ran() * 7.56;
        r[2] = rng.ran() * 28.56;

        // sample omega
        costheta = 1.0 - 2.0 * rng.ran();
        phi      = profugus::constants::two_pi * rng.ran();
        sintheta = sqrt(1.0 - costheta * costheta);

        omega[0] = sintheta * cos(phi);
        omega[1] = sintheta * sin(phi);
        omega[2] = costheta;

        // initialize track
        geometry.initialize(r, omega, state);
        EXPECT_EQ(INSIDE, geometry.boundary_state(state));

        while (geometry.boundary_state(state) == INSIDE)
        {
            // get distance-to-boundary
            d = geometry.distance_to_boundary(state);

            // update position of particle and cross the surface
            geometry.move_to_surface(state);

            if (do_output)
                if (state.d_r[2] < 14.2799 &&
                    state.escaping_face != State::MINUS_Z)
                    csites << state.d_r[0] << "\t" << state.d_r[1] << endl;
        }

        if (state.escaping_face == State::MINUS_X)
            face_bin[0]++;
        else if (state.escaping_face == State::PLUS_X)
            face_bin[1]++;
        else if (state.escaping_face == State::MINUS_Y)
            face_bin[2]++;
        else if (state.escaping_face == State::PLUS_Y)
            face_bin[3]++;
        else if (state.escaping_face == State::MINUS_Z)
            face_bin[4]++;
        else if (state.escaping_face == State::PLUS_Z)
            face_bin[5]++;
    }

    EXPECT_TRUE(face_bin[0] + face_bin[1] + face_bin[2] +
              face_bin[3] + face_bin[4] + face_bin[5] == Np);

    double xyf  = 28.56 * 7.56;
    double zf   = 7.56 * 7.56;
    double area = 4 * xyf + 2 * zf;
    double Npx  = static_cast<double>(Np);
    double lox  = face_bin[0] / Npx;
    double hix  = face_bin[1] / Npx;
    double loy  = face_bin[2] / Npx;
    double hiy  = face_bin[3] / Npx;
    double loz  = face_bin[4] / Npx;
    double hiz  = face_bin[5] / Npx;

    EXPECT_SOFTEQ(lox, xyf / area, 0.01);
    EXPECT_SOFTEQ(hix, xyf / area, 0.03);
    EXPECT_SOFTEQ(loy, xyf / area, 0.03);
    EXPECT_SOFTEQ(hiy, xyf / area, 0.02);
    EXPECT_SOFTEQ(loz, zf / area, 0.04);
    EXPECT_SOFTEQ(hiz, zf / area, 0.02);

    cout.precision(5);
    cout << endl;
    cout << "Low  X leakage = "
         << setw(8) << lox << " ("
         << setw(8) << xyf / area << ")" << endl;
    cout << "High X leakage = "
         << setw(8) << hix << " ("
         << setw(8) << xyf / area << ")" << endl;
    cout << "Low  Y leakage = "
         << setw(8) << loy << " ("
         << setw(8) << xyf / area << ")" << endl;
    cout << "High Y leakage = "
         << setw(8) << hiy << " ("
         << setw(8) << xyf / area << ")" << endl;
    cout << "Low  Z leakage = "
         << setw(8) << loz << " ("
         << setw(8) << zf / area << ")" << endl;
    cout << "High Z leakage = "
         << setw(8) << hiz << " ("
         << setw(8) << zf / area << ")" << endl;
    cout << endl;

    csites.close();

    // Check cell bounding boxes
    {
        using def::I; using def::J; using def::K;
        def::Space_Vector lower, upper;

        // Cell 0
        auto bbox = rtk_core.get_cell_extents(0);
        lower = bbox.lower();
        upper = bbox.upper();

        EXPECT_SOFT_EQ( 0.0,   lower[I] );
        EXPECT_SOFT_EQ( 1.26,  upper[I] );
        EXPECT_SOFT_EQ( 0.0,   lower[J] );
        EXPECT_SOFT_EQ( 1.26,  upper[J] );
        EXPECT_SOFT_EQ( 0.0,   lower[K] );
        EXPECT_SOFT_EQ( 14.28, upper[K] );

        // Cell 12, Core Location (1,0,0), Lattice Location (0, 1)
        bbox = rtk_core.get_cell_extents(12);
        lower = bbox.lower();
        upper = bbox.upper();

        EXPECT_SOFT_EQ( 2.52,  lower[I] );
        EXPECT_SOFT_EQ( 3.78,  upper[I] );
        EXPECT_SOFT_EQ( 1.26,  lower[J] );
        EXPECT_SOFT_EQ( 2.52,  upper[J] );
        EXPECT_SOFT_EQ( 0.0,   lower[K] );
        EXPECT_SOFT_EQ( 14.28, upper[K] );

        // Cell 35, Core Location (1,2,0)
        bbox = rtk_core.get_cell_extents(35);
        lower = bbox.lower();
        upper = bbox.upper();

        EXPECT_SOFT_EQ( 2.52,  lower[I] );
        EXPECT_SOFT_EQ( 5.04,  upper[I] );
        EXPECT_SOFT_EQ( 5.04,  lower[J] );
        EXPECT_SOFT_EQ( 7.56,  upper[J] );
        EXPECT_SOFT_EQ( 0.0,   lower[K] );
        EXPECT_SOFT_EQ( 14.28, upper[K] );

        // Cell 42, Core Location (2,1,1)
        bbox = rtk_core.get_cell_extents(42);
        lower = bbox.lower();
        upper = bbox.upper();

        EXPECT_SOFT_EQ( 5.04,  lower[I] );
        EXPECT_SOFT_EQ( 7.56,  upper[I] );
        EXPECT_SOFT_EQ( 2.52,  lower[J] );
        EXPECT_SOFT_EQ( 5.04,  upper[J] );
        EXPECT_SOFT_EQ( 14.28, lower[K] );
        EXPECT_SOFT_EQ( 28.56, upper[K] );
    }
}

//---------------------------------------------------------------------------//

TEST(Core, Reflecting)
{
    // 2 fuel pin types
    SP_Pin_Cell pin1(make_shared<Pin_Cell_t>(1, 0.54, 3, 1.26, 14.28));
    SP_Pin_Cell pin2(make_shared<Pin_Cell_t>(2, 0.54, 3, 1.26, 14.28));

    // water pin
    SP_Pin_Cell box(make_shared<Pin_Cell_t>(3, 2.52, 14.28));

    // 3 lattices (fuel 1, 2, and water)
    SP_Lattice lat1(make_shared<Lattice_t>(2, 2, 1, 1));
    SP_Lattice lat2(make_shared<Lattice_t>(2, 2, 1, 1));
    SP_Lattice lat3(make_shared<Lattice_t>(1, 1, 1, 1));

    // lattice assignments
    lat1->assign_object(pin1, 0);
    lat2->assign_object(pin2, 0);
    lat3->assign_object(box, 0);

    // complete the lattices
    lat1->complete(0.0, 0.0, 0.0);
    lat2->complete(0.0, 0.0, 0.0);
    lat3->complete(0.0, 0.0, 0.0);

    // make core (3x3x2 with 4 objects, object 0 unassigned)
    SP_Core core(make_shared<Core_t>(3, 3, 2, 4));
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

    // set reflected faces
    vector<int> refl(6, 0);
    refl[0] = 1; refl[2] = 1; refl[4] = 1;
    core->set_reflecting(refl);

    // complete the core
    core->complete(0.0, 0.0, 0.0);

    // check lattice
    EXPECT_TRUE(soft_equiv(core->pitch(X), 7.56));
    EXPECT_TRUE(soft_equiv(core->pitch(Y), 7.56));
    EXPECT_TRUE(soft_equiv(core->height(), 28.56));

    // build the RTK Core
    Core_Geometry rtk_core(core);
    Geometry &geometry = rtk_core;

    // make a random number generator
    profugus::RNG_Control control(seed);
    auto rng = control.rng();

    // geometry variables
    double costheta, sintheta, phi;
    Vector r, omega;
    State  state;
    double d;
    int    Np = 10000;

    int face_bin[6] = {0};
    int refl_bin[6] = {0};

    // sample Np tracks
    for (int n = 0; n < Np; ++n)
    {
        // sample x,y,z randomly
        r[0] = rng.ran() * 7.56;
        r[1] = rng.ran() * 7.56;
        r[2] = rng.ran() * 28.56;

        // sample omega
        costheta = 1.0 - 2.0 * rng.ran();
        phi      = profugus::constants::two_pi * rng.ran();
        sintheta = sqrt(1.0 - costheta * costheta);

        omega[0] = sintheta * cos(phi);
        omega[1] = sintheta * sin(phi);
        omega[2] = costheta;

        // initialize track
        geometry.initialize(r, omega, state);
        EXPECT_EQ(INSIDE, geometry.boundary_state(state));

        // continue flag
        bool done = false;

        while (!done)
        {
            // get distance-to-boundary
            d = geometry.distance_to_boundary(state);
            EXPECT_TRUE(geometry.boundary_state(state) != REFLECT);

            // update position of particle to the surface and process it through
            geometry.move_to_surface(state);

            // if the particle is reflected then do the reflection
            if (geometry.boundary_state(state) == REFLECT)
            {
                if (state.exiting_face == State::MINUS_X)
                    refl_bin[0]++;
                else if (state.exiting_face == State::PLUS_X)
                    refl_bin[1]++;
                else if (state.exiting_face == State::MINUS_Y)
                    refl_bin[2]++;
                else if (state.exiting_face == State::PLUS_Y)
                    refl_bin[3]++;
                else if (state.exiting_face == State::MINUS_Z)
                    refl_bin[4]++;
                else if (state.exiting_face == State::PLUS_Z)
                    refl_bin[5]++;

                // reflect the particle
                EXPECT_TRUE(geometry.reflect(state));
            }

            // terminate on escape
            if (geometry.boundary_state(state) == OUTSIDE)
            {
                done = true;
            }
        }

        if (state.escaping_face == State::MINUS_X)
            face_bin[0]++;
        else if (state.escaping_face == State::PLUS_X)
            face_bin[1]++;
        else if (state.escaping_face == State::MINUS_Y)
            face_bin[2]++;
        else if (state.escaping_face == State::PLUS_Y)
            face_bin[3]++;
        else if (state.escaping_face == State::MINUS_Z)
            face_bin[4]++;
        else if (state.escaping_face == State::PLUS_Z)
            face_bin[5]++;
    }

    EXPECT_TRUE(face_bin[0] + face_bin[1] + face_bin[2] +
              face_bin[3] + face_bin[4] + face_bin[5] == Np);

    EXPECT_EQ(0, face_bin[0]);
    EXPECT_EQ(0, refl_bin[1]);
    EXPECT_EQ(0, face_bin[2]);
    EXPECT_EQ(0, refl_bin[3]);
    EXPECT_EQ(0, face_bin[4]);
    EXPECT_EQ(0, refl_bin[5]);

    // heuristicly stored data
    EXPECT_EQ(2992, refl_bin[0]);
    EXPECT_EQ(4563, face_bin[1]);
    EXPECT_EQ(2989, refl_bin[2]);
    EXPECT_EQ(4281, face_bin[3]);
    EXPECT_EQ(1096, refl_bin[4]);
    EXPECT_EQ(1156, face_bin[5]);

    cout.precision(5);
    cout << endl;
    cout << "Low  X leakage    = " << setw(8) << face_bin[0] << endl;
    cout << "High X leakage    = " << setw(8) << face_bin[1] << endl;
    cout << "Low  Y leakage    = " << setw(8) << face_bin[2] << endl;
    cout << "High Y leakage    = " << setw(8) << face_bin[3] << endl;
    cout << "Low  Z leakage    = " << setw(8) << face_bin[4] << endl;
    cout << "High Z leakage    = " << setw(8) << face_bin[5] << endl;
    cout << endl;

    cout << "Low  X reflection = " << setw(8) << refl_bin[0] << endl;
    cout << "High X reflection = " << setw(8) << refl_bin[1] << endl;
    cout << "Low  Y reflection = " << setw(8) << refl_bin[2] << endl;
    cout << "High Y reflection = " << setw(8) << refl_bin[3] << endl;
    cout << "Low  Z reflection = " << setw(8) << refl_bin[4] << endl;
    cout << "High Z reflection = " << setw(8) << refl_bin[5] << endl;
    cout << endl;
}

//---------------------------------------------------------------------------//
// See support/bwr.png for core figure showing particle path.

TEST(Bwr, Lattice)
{
    // make bwr lattice
    SP_Core bwr;
    {
        // pin shells
        vector<double> r(2, 0.0);
        vector<int>    rid(2, 1);
        r[0] = 0.125;
        r[1] = 0.25;

        // cells
        SP_Pin_Cell pin(make_shared<Pin_Cell_t>(rid, r, 5, 0.75, 14.00, 4));
        SP_Pin_Cell plug(make_shared<Pin_Cell_t>(2, 0.6, 5, 1.5, 14.00, 4));

        // gaps
        SP_Pin_Cell g0(make_shared<Pin_Cell_t>(5, 0.25, 14.00));      // corner
        SP_Pin_Cell g1(make_shared<Pin_Cell_t>(5, 0.25, 1.5, 14.00)); // x-edge
        SP_Pin_Cell g2(make_shared<Pin_Cell_t>(5, 1.5, 0.25, 14.00)); // y-edge

        // 2x2 pin lattice
        SP_Lattice lat1(make_shared<Lattice_t>(2, 2, 1, 1));
        lat1->assign_object(pin, 0);
        lat1->complete(0.0, 0.0, 0.0);
        EXPECT_EQ(48, lat1->num_cells());

        // plug lattice
        SP_Lattice lat2(make_shared<Lattice_t>(1, 1, 1, 1));
        lat2->assign_object(plug, 0);
        lat2->complete(0.0, 0.0, 0.0);
        EXPECT_EQ(8, lat2->num_cells());

        // corner gap lattice
        SP_Lattice lat3(make_shared<Lattice_t>(1, 1, 1, 1));
        lat3->assign_object(g0, 0);
        lat3->complete(0.0, 0.0, 0.0);
        EXPECT_EQ(1, lat3->num_cells());

        // x-edge gap lattice
        SP_Lattice lat4(make_shared<Lattice_t>(1, 1, 1, 1));
        lat4->assign_object(g1, 0);
        lat4->complete(0.0, 0.0, 0.0);
        EXPECT_EQ(1, lat4->num_cells());

        // y-edge gap lattice
        SP_Lattice lat5(make_shared<Lattice_t>(1, 1, 1, 1));
        lat5->assign_object(g2, 0);
        lat5->complete(0.0, 0.0, 0.0);
        EXPECT_EQ(1, lat5->num_cells());

        bwr = make_shared<Core_t>(4, 4, 1, 5);

        bwr->assign_object(lat1, 0); // 2x2 pin lattice
        bwr->assign_object(lat2, 1); // water plug
        bwr->assign_object(lat3, 2); // corner gap
        bwr->assign_object(lat4, 3); // x-edge gap
        bwr->assign_object(lat5, 4); // y_edge gap

        bwr->id(0, 0, 0) = 2;
        bwr->id(1, 0, 0) = 4;
        bwr->id(2, 0, 0) = 4;
        bwr->id(3, 0, 0) = 2;

        bwr->id(0, 1, 0) = 3;
        bwr->id(1, 1, 0) = 1;
        bwr->id(2, 1, 0) = 0;
        bwr->id(3, 1, 0) = 3;

        bwr->id(0, 2, 0) = 3;
        bwr->id(1, 2, 0) = 0;
        bwr->id(2, 2, 0) = 1;
        bwr->id(3, 2, 0) = 3;

        bwr->id(0, 3, 0) = 2;
        bwr->id(1, 3, 0) = 4;
        bwr->id(2, 3, 0) = 4;
        bwr->id(3, 3, 0) = 2;

        // complete the bwr lattice
        bwr->complete(0.0, 0.0, 0.0);

        EXPECT_TRUE(soft_equiv(bwr->pitch(X), 3.5));
        EXPECT_TRUE(soft_equiv(bwr->pitch(Y), 3.5));
    }

    // build the RTK Bwr
    Core_Geometry bg(bwr);
    Geometry &geometry = bg;

    cout << endl;

    Vector r, omega;
    State  state;
    double eps = 1.0e-6, d;

    // initial point
    r[0] = 1.3;
    r[1] = 3.5;
    r[2] = 0.0;

    // direction
    omega[0] = 2.2;
    omega[1] = -3.2;
    omega[2] = 0.1;
    profugus::vector_normalize(omega);

    // initialize the geoemtry
    bg.initialize(r, omega, state);

    int ref[] = {121, 101, 107, 114, 113, 117, 118, 116,
                 54, 53, 52, 58, 59, 56, 57, 27, 61};
    int k = 0;

    cout.precision(5);
    while (bg.boundary_state(state) != OUTSIDE)
    {
        d = bg.distance_to_boundary(state);

        int cell = bg.cell(state);
        int mat  = bg.matid(state);

        EXPECT_EQ(ref[k], cell);
        k++;

        bg.move_to_surface(state);

        cout << "Transporting " << fixed << setw(10) << d
             << " in cell " << setw(4) << cell
             << " through material " << mat
             << " to " << setw(10) << bg.position(state)[0]
             << setw(10) << bg.position(state)[1]
             << setw(10) << bg.position(state)[2] << endl;
    }

    cout << endl;
}

//---------------------------------------------------------------------------//
// See support/bwr.png for core figure showing particle path.

TEST(Bwr, Symmetry)
{
    // make bwr lattice
    SP_Core bwr;
    {
        // pin shells
        vector<double> r(2, 0.0);
        vector<int>    rid(2, 1);
        r[0] = 0.125;
        r[1] = 0.25;

        // cells
        SP_Pin_Cell pin(make_shared<Pin_Cell_t>(rid, r, 5, 0.75, 14.00, 4));
        SP_Pin_Cell plug(make_shared<Pin_Cell_t>(2, 0.6, 5, 1.5, 14.00, 4));

        // gaps
        SP_Pin_Cell g0(make_shared<Pin_Cell_t>(5, 0.25, 14.00));      // corner
        SP_Pin_Cell g1(make_shared<Pin_Cell_t>(5, 0.25, 1.5, 14.00)); // x-edge
        SP_Pin_Cell g2(make_shared<Pin_Cell_t>(5, 1.5, 0.25, 14.00)); // y-edge

        // 2x2 pin lattice
        SP_Lattice lat1(make_shared<Lattice_t>(2, 2, 1, 1));
        lat1->assign_object(pin, 0);
        lat1->complete(0.0, 0.0, 0.0);
        EXPECT_EQ(48, lat1->num_cells());

        // plug lattice
        SP_Lattice lat2(make_shared<Lattice_t>(1, 1, 1, 1));
        lat2->assign_object(plug, 0);
        lat2->complete(0.0, 0.0, 0.0);
        EXPECT_EQ(8, lat2->num_cells());

        // corner gap lattice
        SP_Lattice lat3(make_shared<Lattice_t>(1, 1, 1, 1));
        lat3->assign_object(g0, 0);
        lat3->complete(0.0, 0.0, 0.0);
        EXPECT_EQ(1, lat3->num_cells());

        // x-edge gap lattice
        SP_Lattice lat4(make_shared<Lattice_t>(1, 1, 1, 1));
        lat4->assign_object(g1, 0);
        lat4->complete(0.0, 0.0, 0.0);
        EXPECT_EQ(1, lat4->num_cells());

        // y-edge gap lattice
        SP_Lattice lat5(make_shared<Lattice_t>(1, 1, 1, 1));
        lat5->assign_object(g2, 0);
        lat5->complete(0.0, 0.0, 0.0);
        EXPECT_EQ(1, lat5->num_cells());

        bwr = make_shared<Core_t>(4, 4, 1, 5);

        bwr->assign_object(lat1, 0); // 2x2 pin lattice
        bwr->assign_object(lat2, 1); // water plug
        bwr->assign_object(lat3, 2); // corner gap
        bwr->assign_object(lat4, 3); // x-edge gap
        bwr->assign_object(lat5, 4); // y_edge gap

        bwr->id(0, 0, 0) = 2;
        bwr->id(1, 0, 0) = 4;
        bwr->id(2, 0, 0) = 4;
        bwr->id(3, 0, 0) = 2;

        bwr->id(0, 1, 0) = 3;
        bwr->id(1, 1, 0) = 1;
        bwr->id(2, 1, 0) = 0;
        bwr->id(3, 1, 0) = 3;

        bwr->id(0, 2, 0) = 3;
        bwr->id(1, 2, 0) = 0;
        bwr->id(2, 2, 0) = 1;
        bwr->id(3, 2, 0) = 3;

        bwr->id(0, 3, 0) = 2;
        bwr->id(1, 3, 0) = 4;
        bwr->id(2, 3, 0) = 4;
        bwr->id(3, 3, 0) = 2;

        // complete the bwr lattice
        bwr->complete(0.0, 0.0, 0.0);

        EXPECT_TRUE(soft_equiv(bwr->pitch(X), 3.5));
        EXPECT_TRUE(soft_equiv(bwr->pitch(Y), 3.5));
    }
}

//---------------------------------------------------------------------------//
//                        end of tstRTK_Core.cc
//---------------------------------------------------------------------------//
