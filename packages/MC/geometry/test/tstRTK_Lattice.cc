//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/test/tstRTK_Lattice.cc
 * \author Thomas M. Evans
 * \date   Wed Dec 22 13:32:47 2010
 * \brief  RTK_Lattice unit test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/nemesis_gtest.hh"
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include <geometry/config.h>

#include "utils/Definitions.hh"
#include "utils/Constants.hh"
#include "utils/Vector_Functions.hh"
#include "geometry/Definitions.hh"
#include "../RTK_Geometry.hh"

#ifdef USE_MC
#include "mc/RNG_Control.hh"
#endif

using namespace std;

using def::X;
using def::Y;
using def::Z;
using nemesis::constants::pi;

typedef denovo::RTK_Lattice         Lattice_Geometry;
typedef Lattice_Geometry            Geometry;
typedef Lattice_Geometry::Array_t   Lattice_t;
typedef Lattice_t::Object_t         Pin_Cell_t;
typedef Lattice_Geometry::SP_Array  SP_Lattice;
typedef Lattice_t::SP_Object        SP_Pin_Cell;

typedef Geometry::Space_Vector Vector;
typedef Geometry::Geo_State_t  State;

using denovo::geometry::OUTSIDE;
using denovo::geometry::INSIDE;


int seed = 4305834;


bool do_output = false;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Single, ShellLattice)
{
    // Fuel region = 1
    // Fuel region = 2
    // Moderator   = 5

    // 3x3x1 lattice, 2 fuel pin types, 1 box
    SP_Pin_Cell pin1(new Pin_Cell_t(1, 0.54, 5, 1.26, 14.28));
    SP_Pin_Cell pin2(new Pin_Cell_t(2, 0.54, 5, 1.26, 14.28));
    SP_Pin_Cell box(new Pin_Cell_t(5, 1.26, 14.28));

    SP_Lattice lat(new Lattice_t(3, 3, 1, 3));

    // assign pins
    lat->assign_object(pin1, 0);
    lat->assign_object(pin2, 1);
    lat->assign_object(box,  2);

    // arrange pin-cells in lattice
    lat->id(0, 0, 0) = 0; // pin 1
    lat->id(1, 0, 0) = 0; // pin 1
    lat->id(2, 0, 0) = 0; // pin 1
    lat->id(0, 1, 0) = 0; // pin 1
    lat->id(1, 1, 0) = 2; // box
    lat->id(2, 1, 0) = 1; // pin 2
    lat->id(0, 2, 0) = 0; // pin 1
    lat->id(1, 2, 0) = 1; // pin 2
    lat->id(2, 2, 0) = 1; // pin 2

    // complete lattice
    lat->complete(0.0, 0.0, 0.0);

    // build the lattice geometry
    Lattice_Geometry lattice(lat);

    // test the geometry interface
    State  state;
    double d;
    int    m, step = 0;
    Geometry &geometry = lattice;

    Vector pos, dir;

    cout.precision(3);
    {
        double dref[] = {0.009308327, 0.092021331, 0.091669426, 1.093497688,
                         0.091586583};

        // start a track
        pos = Vector(  1.36,   0.60,   3.80);
        dir = Vector( -0.986877895121,   0.039161821235,   0.156647284940);
        geometry.initialize(pos, dir, state);
        EXPECT_EQ(1, geometry.matid(state));
        EXPECT_EQ(INSIDE, geometry.boundary_state(state));

        cout << setw(3) << step << ": starting at (" << fixed
             << setw(6) << state.d_r[0] << ","
             << setw(6) << state.d_r[1] << ","
             << setw(6) << state.d_r[2]
             << ") in direction ("
             << setw(6) << state.d_dir[0] << ","
             << setw(6) << state.d_dir[1] << ","
             << setw(6) << state.d_dir[2]
             << ")" << endl;

        while (geometry.boundary_state(state) == INSIDE)
        {
            d = geometry.distance_to_boundary(state);
            m = geometry.matid(state);

            EXPECT_SOFTEQ(d, dref[step], 1.0e-6);

            // move to the next surface
            geometry.move_to_surface(state);

            // move particle to next surface
            cout << setw(3) << ++step << ": moving particle d = "
                 << fixed << setw(7)
                 << d << " through mat " << m << " to ("
                 << setw(8) << state.d_r[0] << ","
                 << setw(8) << state.d_r[1] << ","
                 << setw(8) << state.d_r[2] << ")" << endl;
        }
        cout << endl;
        EXPECT_EQ(State::MINUS_X, state.escaping_face);
    }

    step = 0;

    {
        double dref[] = {0.573902948, 0.143274851, 0.645971956, 1.318471580,
                         0.382091481, 0.206180701, 1.073886420, 0.302284933};

        // start a track
        pos = Vector(  0.25,   0.80,   3.10);
        dir = Vector(  0.740931065010,   0.641403011501,   0.199056107018);
        geometry.initialize(pos, dir, state);
        EXPECT_EQ(1, geometry.matid(state));
        EXPECT_EQ(INSIDE, geometry.boundary_state(state));

        EXPECT_EQ(0.25, geometry.position(state)[0]);
        EXPECT_EQ(0.80, geometry.position(state)[1]);
        EXPECT_EQ(3.10, geometry.position(state)[2]);

        cout << setw(3) << step << ": starting at (" << fixed
             << setw(6) << state.d_r[0] << ","
             << setw(6) << state.d_r[1] << ","
             << setw(6) << state.d_r[2]
             << ") in direction ("
             << setw(6) << state.d_dir[0] << ","
             << setw(6) << state.d_dir[1] << ","
             << setw(6) << state.d_dir[2]
             << ")" << endl;

        while (geometry.boundary_state(state) == INSIDE)
        {
            d = geometry.distance_to_boundary(state);
            m = geometry.matid(state);

            EXPECT_SOFTEQ(d, dref[step], 1.0e-6);

            // move to next surface
            geometry.move_to_surface(state);

            // move particle to next surface
            cout << setw(3) << ++step << ": moving particle d = "
                 << fixed << setw(7)
                 << d << " through mat " << m << " to ("
                 << setw(8) << state.d_r[0] << ","
                 << setw(8) << state.d_r[1] << ","
                 << setw(8) << state.d_r[2] << ")" << endl;
        }
        cout << endl;
        EXPECT_EQ(State::PLUS_Y, state.escaping_face);
    }

    // change direction
    double v = 0.57735026918962584;
    geometry.change_direction(Vector(v, v, v), state);
    EXPECT_SOFTEQ(state.d_dir[0], v, 1.e-12);
    EXPECT_SOFTEQ(state.d_dir[1], v, 1.e-12);
    EXPECT_SOFTEQ(state.d_dir[2], v, 1.e-12);

    // reflect off of various faces (we aren't moving explicitly to the faces,
    // but it doesn't check position so we get away with it for the purposes
    // of this test)
    geometry.change_direction(Vector(1.0, 0.0, 0.0), state);
    state.exiting_face = State::PLUS_X;
    EXPECT_TRUE(geometry.reflect(state));
    EXPECT_SOFTEQ(state.d_dir[0], -1.0, 1.e-12);
    EXPECT_SOFTEQ(state.d_dir[1], 0.0, 1.e-12);
    EXPECT_SOFTEQ(state.d_dir[2], 0.0, 1.e-12);

    geometry.change_direction(Vector(-1.0, 1.0, -2.0), state);
    double a = state.d_dir[X], b = state.d_dir[Y], c = state.d_dir[Z];
    state.exiting_face = State::MINUS_X;
    geometry.reflect(state);
    EXPECT_SOFTEQ(state.d_dir[0], -a, 1.e-12);
    EXPECT_SOFTEQ(state.d_dir[1], b, 1.e-12);
    EXPECT_SOFTEQ(state.d_dir[2], c, 1.e-12);

    double n   = 1.0/sqrt(3.0);
    double dot = n * (-a + b + c);
    EXPECT_SOFTEQ(state.d_dir[0], -a - 2.0 * n * dot, 1.e-12);
    EXPECT_SOFTEQ(state.d_dir[1], b - 2.0 * n * dot, 1.e-12);
    EXPECT_SOFTEQ(state.d_dir[2], c - 2.0 * n * dot, 1.e-12);

    state.exiting_face = State::INTERNAL;
    EXPECT_TRUE(!geometry.reflect(state));
    EXPECT_SOFTEQ(state.d_dir[0], -a - 2.0 * n * dot, 1.e-12);
    EXPECT_SOFTEQ(state.d_dir[1], b - 2.0 * n * dot, 1.e-12);
    EXPECT_SOFTEQ(state.d_dir[2], c - 2.0 * n * dot, 1.e-12);

    // test change-direction through (theta, phi)
    double costheta = cos(2.0/3.0);
    double phi      = 4.0 * pi / 3.0;
    geometry.change_direction(Vector(1.0, 1.0, 1.0), state);
    geometry.change_direction(costheta, phi, state);
    EXPECT_SOFTEQ(state.d_dir[0], 0.70618063654, 1.e-12);
    EXPECT_SOFTEQ(state.d_dir[1], -0.0511646083931, 1.e-12);
    EXPECT_SOFTEQ(state.d_dir[2], 0.70618063654, 1.e-12);

    costheta = cos(-1.0/3.0);
    phi      = 4.5 * pi / 3.0;
    geometry.change_direction(costheta, phi, state);
    {
        Vector omega = geometry.direction(state);
        EXPECT_SOFTEQ(omega[0], 0.643666175435, 1.0e-11);
        EXPECT_SOFTEQ(omega[1], -0.374687631211, 1.0e-11);
        EXPECT_SOFTEQ(omega[2], 0.667310297851, 1.0e-11);
    }

    // ut.passes("Lattice of single-shell pins correctly represented.");
}

//---------------------------------------------------------------------------//

TEST(Lattice, Heuristic)
{
#ifdef USE_MC
    // Fuel region = 1
    // Fuel region = 2
    // Moderator   = 5

    // 3x3x1 lattice, 2 fuel pin types, 1 box
    SP_Pin_Cell pin1(new Pin_Cell_t(1, 0.54, 5, 1.26, 14.28));
    SP_Pin_Cell pin2(new Pin_Cell_t(2, 0.54, 5, 1.26, 14.28));
    SP_Pin_Cell box(new Pin_Cell_t(5, 1.26, 14.28));

    SP_Lattice lat(new Lattice_t(3, 3, 1, 3));

    // assign pins
    lat->assign_object(pin1, 0);
    lat->assign_object(pin2, 1);
    lat->assign_object(box,  2);

    // arrange pin-cells in lattice
    lat->id(0, 0, 0) = 0; // pin 1
    lat->id(1, 0, 0) = 0; // pin 1
    lat->id(2, 0, 0) = 0; // pin 1
    lat->id(0, 1, 0) = 0; // pin 1
    lat->id(1, 1, 0) = 2; // box
    lat->id(2, 1, 0) = 1; // pin 2
    lat->id(0, 2, 0) = 0; // pin 1
    lat->id(1, 2, 0) = 1; // pin 2
    lat->id(2, 2, 0) = 1; // pin 2

    // complete lattice
    lat->complete(0.0, 0.0, 0.0);

    // build the lattice geometry
    Lattice_Geometry lattice(lat);
    Geometry &geometry = lattice;

    // make a random number generator
    mc::RNG_Control control(seed);
    mc::RNG_Control::RNG rng = control.rng();

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
        // sample x,y,z
        r[0] = rng.ran() * 3.78;
        r[1] = rng.ran() * 3.78;
        r[2] = rng.ran() * 14.28;

        // sample omega
        costheta = 1.0 - 2.0 * rng.ran();
        phi      = nemesis::constants::two_pi * rng.ran();
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
            EXPECT_SOFTEQ(d, state.dist_to_next_region, 1.e-12);

            // update position of particle and cross surface
            geometry.move_to_surface(state);

            if (do_output)
                if (state.escaping_face != State::PLUS_Z &&
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

    double xyf  = 14.28 * 3.78;
    double zf   = 3.78 * 3.78;
    double area = 4 * xyf + 2 * zf;
    double Npx  = static_cast<double>(Np);
    double lox  = face_bin[0] / Npx;
    double hix  = face_bin[1] / Npx;
    double loy  = face_bin[2] / Npx;
    double hiy  = face_bin[3] / Npx;
    double loz  = face_bin[4] / Npx;
    double hiz  = face_bin[5] / Npx;

    EXPECT_SOFTEQ(lox, xyf / area, 0.01);
    EXPECT_SOFTEQ(hix, xyf / area, 0.05);
    EXPECT_SOFTEQ(loy, xyf / area, 0.05);
    EXPECT_SOFTEQ(hiy, xyf / area, 0.05);
    EXPECT_SOFTEQ(loz, zf / area, 0.05);
    EXPECT_SOFTEQ(hiz, zf / area, 0.05);

    cout.precision(5);
    cout << fixed;
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

    // ut.passes("Finished heuristic tracking through lattice.");
#else
    // ut.passes("Heuristic tests need mc package for RNG.");
#endif
}

//---------------------------------------------------------------------------//

TEST(Multisegment, Pin)
{
#ifdef USE_MC

    // make pin with clad and 4 segments
    vector<int>    ids(3, 0);
    vector<double> rad(3, 0);
    ids[0] = 1;
    rad[0] = 0.27;
    ids[1] = 2;
    rad[1] = 0.486;
    ids[2] = 5;
    rad[2] = 0.54;

    // 1x1x1 lattice, 1 fuel pin with 3 shells, 4 regions, 4 segments
    SP_Pin_Cell pin(new Pin_Cell_t(ids, rad, 10, 1.26, 14.28, 4));

    SP_Lattice lat(new Lattice_t(1, 1, 1, 1));

    // assign pins
    lat->assign_object(pin, 0);

    // complete lattice
    lat->complete(0.0, 0.0, 0.0);

    // build the lattice geometry
    Lattice_Geometry lattice(lat);
    Geometry &geometry = lattice;

    // make a random number generator
    mc::RNG_Control control(seed);
    mc::RNG_Control::RNG rng = control.rng();

    // geometry variables
    double costheta, sintheta, phi;
    Vector r, omega;
    State  state;
    double dbnd, dcol, d;
    int    Np = 20000;
    int    mid;

    double xs[] = {0.0, 0.5, 0.6, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.7};

    int face_bin[6] = {0};

    cout.precision(6);

    // volume tallies
    double v[16] = {0.0};
    double vt    = 0.0;

    // reference volumes
    double Vr[4] = {0.0};

    // sample Np tracks
    {
        for (int n = 0; n < Np; ++n)
        {
            // sample y,z on low x face
            r[0] = 1.0e-12;
            r[1] = rng.ran() * 1.26;
            r[2] = rng.ran() * 14.28;

            omega[0] = 1.0;
            omega[1] = 0.0;
            omega[2] = 0.0;

            // initialize track
            geometry.initialize(r, omega, state);
            EXPECT_EQ(INSIDE, geometry.boundary_state(state));

            while (geometry.boundary_state(state) == INSIDE)
            {
                // get distance-to-boundary
                dbnd = geometry.distance_to_boundary(state);
                EXPECT_SOFTEQ(dbnd, state.dist_to_next_region, 1.e-12);

                // get distance to collision
                mid  = geometry.matid(state);
                dcol = -log(rng.ran()) / xs[mid];
                EXPECT_TRUE(xs[mid] != 0.0);

                EXPECT_TRUE(geometry.cell(state) ==
                          state.region + state.segment * 4);
                EXPECT_TRUE(geometry.cell(state) < 16);

                if (dcol < dbnd)
                {
                    // tally the total path
                    vt += dcol;

                    // tally the path through this particular cell
                    v[geometry.cell(state)] += dcol;

                    // move to the point
                    geometry.move_to_point(dcol, state);
                }
                else
                {
                    // tally the total path
                    vt += dbnd;

                    // tally the path through this particular cell
                    v[geometry.cell(state)] += dbnd;

                    // update position of particle and cross surface
                    geometry.move_to_surface(state);
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

        double xyf  = 14.28 * 1.26;
        double zf   = 1.26 * 1.26;
        double area = 4 * xyf + 2 * zf;
        double Npx  = static_cast<double>(Np);
        double lox  = face_bin[0] / Npx;
        double hix  = face_bin[1] / Npx;
        double loy  = face_bin[2] / Npx;
        double hiy  = face_bin[3] / Npx;
        double loz  = face_bin[4] / Npx;
        double hiz  = face_bin[5] / Npx;

        EXPECT_SOFTEQ(hix, 1.0, 1.e-12);
        EXPECT_EQ(0.0, lox);
        EXPECT_EQ(0.0, loy);
        EXPECT_EQ(0.0, loz);
        EXPECT_EQ(0.0, hiy);
        EXPECT_EQ(0.0, hiz);

        cout.precision(5);
        cout << fixed;
        cout << endl;

        // region 0
        Vr[0] = pi*rad[0]*rad[0] * 14.28 * 0.25;

        // region 1
        Vr[1] = pi*(rad[1]*rad[1] - rad[0]*rad[0]) * 14.28 * 0.25;

        // region 2
        Vr[2] = pi*(rad[2]*rad[2] - rad[1]*rad[1]) * 14.28 * 0.25;

        // region 3
        Vr[3] = (1.26*1.26 - pi*rad[2]*rad[2]) * 14.28 * 0.25;

        // volume tallies
        double Vtot = 1.26 * 1.26 * 14.28;
        double sum  = 0.0;
        for (int s = 0; s < 4; ++s)
        {
            for (int r = 0; r < 4; ++r)
            {
                int cell = r + s * 4;

                double Vol = v[cell] / vt * Vtot;
                double err = fabs(Vol - Vr[r]) / Vr[r];
                sum       += Vr[r];

                cout << "Volume in cell " << setw(4) << cell
                     << "(" << r << "/" << s << ") ="
                     << fixed << setw(10) << Vol
                     << fixed << setw(10) << Vr[r]
                     << scientific << setw(16) << err << endl;

                EXPECT_TRUE(err < 3.0e-2);
            }
        }
        cout << endl;

        EXPECT_SOFTEQ(sum, Vtot, 1.e-12);
    }

    // check non-orthogonal angle
    {
        vt = 0.0;
        fill(v, v + 16, 0.0);

        for (int n = 0; n < Np; ++n)
        {
            // sample face
            if (rng.ran() < 0.5)
            {
                // low x face
                r[0] = 1.0e-12;
                r[1] = rng.ran() * 1.26;
            }
            else
            {
                // low y face
                r[1] = 1.0e-12;
                r[0] = rng.ran() * 1.26;

            }
            r[2] = rng.ran() * 14.28;

            omega[0] = 1.0;
            omega[1] = 1.0;
            omega[2] = 0.0;
            nemesis::vector_normalize(omega);

            // initialize track
            geometry.initialize(r, omega, state);
            EXPECT_EQ(INSIDE, geometry.boundary_state(state));

            while (geometry.boundary_state(state) == INSIDE)
            {
                // get distance-to-boundary
                dbnd = geometry.distance_to_boundary(state);
                EXPECT_SOFTEQ(dbnd, state.dist_to_next_region, 1.e-12);

                // get distance to collision
                mid  = geometry.matid(state);
                dcol = -log(rng.ran()) / xs[mid];
                EXPECT_TRUE(xs[mid] != 0.0);

                EXPECT_TRUE(geometry.cell(state) ==
                          state.region + state.segment * 4);
                EXPECT_TRUE(geometry.cell(state) < 16);

                if (dcol < dbnd)
                {
                    // tally the total path
                    vt += dcol;

                    // tally the path through this particular cell
                    v[geometry.cell(state)] += dcol;

                    // move to the point
                    geometry.move_to_point(dcol, state);
                }
                else
                {
                    // tally the total path
                    vt += dbnd;

                    // tally the path through this particular cell
                    v[geometry.cell(state)] += dbnd;

                    // update position of particle and cross surface
                    geometry.move_to_surface(state);
                }
            }
        }

        cout.precision(5);
        cout << fixed;
        cout << endl;

        // volume tallies
        double Vtot = 1.26 * 1.26 * 14.28;
        double sum  = 0.0;
        for (int s = 0; s < 4; ++s)
        {
            for (int r = 0; r < 4; ++r)
            {
                int cell = r + s * 4;

                double Vol = v[cell] / vt * Vtot;
                double err = fabs(Vol - Vr[r]) / Vr[r];
                sum       += Vr[r];

                cout << "Volume in cell " << setw(4) << cell
                     << "(" << r << "/" << s << ") ="
                     << fixed << setw(10) << Vol
                     << fixed << setw(10) << Vr[r]
                     << scientific << setw(16) << err << endl;

                EXPECT_TRUE(err < 3.0e-2);
            }
        }
        cout << endl;

        EXPECT_SOFTEQ(sum, Vtot, 1.e-12);
    }

    // ut.passes("Finished heuristic tracking through multi-segment pin-cell.");
#else
    // ut.passes("Heuristic tests need mc package for RNG.");
#endif
}

//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//
//                        end of tstRTK_Lattice.cc
//---------------------------------------------------------------------------//
