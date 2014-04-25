//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/tstRNG.cc
 * \author Thomas M. Evans
 * \date   Fri Apr 25 15:05:58 2014
 * \brief  RNG test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "../RNG.hh"

// Pull in C interface
#include "rng/sprng/sprng.h"

using profugus::RNG;
using namespace std;

int seed = 493875348;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

double lfg_ref[] = {
    0.358227372451, 0.814400762544, 0.877915224103, 0.013705751228,
    0.401882104810, 0.276955195427, 0.792529260078, 0.987080323020,
    0.352553167391, 0.859948774993, 0.603194308422, 0.240573628072,
    0.389723295841, 0.981569511369, 0.409877327483, 0.039204148105,
    0.189383281469, 0.865073916312, 0.304706611691, 0.232801632573,
    0.794369955697, 0.475270265306, 0.484947205790, 0.270386943521,
    0.729128285948, 0.823258029158, 0.963868221039, 0.586210907392,
    0.283874949773, 0.626683750055, 0.682020350364, 0.191056284682,
    0.848765828742, 0.389899491974, 0.305126966944, 0.960845022653,
    0.597414355975, 0.853748950574, 0.225555172754, 0.733257881909,
    0.622582092083, 0.945228188222, 0.107676679222, 0.045038984284,
    0.894956223142, 0.272204785684, 0.295165929398, 0.783339838921,
    0.380358306623, 0.075074055619, 0.114718778590, 0.601271568731,
    0.214596651344, 0.164118849436, 0.630805887378, 0.545005429186,
    0.055243134502, 0.228240561446, 0.551247496056, 0.291540662953,
    0.675499341066, 0.705321802045, 0.386783999752, 0.390331348605,
    0.993024707867, 0.500758546674, 0.713764661966, 0.362301054948,
    0.869157130342, 0.415543548917, 0.804086390529, 0.135070412602,
    0.417567174889, 0.962812618860, 0.188322432146, 0.740913905587,
    0.237594532163, 0.540964664200, 0.461585505553, 0.000743970526,
    0.369983235566, 0.375259839368, 0.214822128501, 0.089939152244,
    0.662357829376, 0.997399917348, 0.003135562211, 0.364867607853,
    0.170392839547, 0.666014827604, 0.204853306253, 0.956194240405,
    0.448814065555, 0.130463086639, 0.073960599598, 0.834785113690,
    0.099415951300, 0.014827871356, 0.667904945291, 0.197449235856};

TEST(RNG, lfg)
{
    // data for random numbers
    int *id;
    int num = 5;

    // get two random number seeds
    id = init_sprng(4, num, seed, 1);

    // make a sprng object
    RNG ran(id, 4);
    EXPECT_EQ(4, ran.get_num());
    ran.print();

    for (int i = 0; i < 100; ++i)
    {
        EXPECT_SOFTEQ(lfg_ref[i], ran.ran(), 1.0e-10);
    }
}

TEST(RNG, lfg_uniform_dbl)
{
    // data for random numbers
    int *id;
    int num = 5;

    // get two random number seeds
    id = init_sprng(4, num, seed, 1);

    // make a sprng object
    RNG ran(id, 4);
    EXPECT_EQ(4, ran.get_num());

    for (int i = 0; i < 100; ++i)
    {
        EXPECT_SOFTEQ(lfg_ref[i], ran.uniform<double>(), 1.0e-10);
    }
}

TEST(RNG, lfg_uniform_flt)
{
    // data for random numbers
    int *id;
    int num = 5;

    // get two random number seeds
    id = init_sprng(4, num, seed, 1);

    // make a sprng object
    RNG ran(id, 4);
    EXPECT_EQ(4, ran.get_num());

    for (int i = 0; i < 100; ++i)
    {
        EXPECT_SOFTEQ(lfg_ref[i], ran.uniform<float>(), 1.0e-6);
    }
}
//---------------------------------------------------------------------------//

TEST(RNG, Randomness)
{
    // data for random numbers
    int *id1, *id2, *id3;
    int num = 5;

    // get two random number seeds
    id1 = init_sprng(0, num, seed, 1);
    id2 = init_sprng(1, num, seed, 1);
    id3 = init_sprng(0, num, seed, 1);

    RNG ran1(id1, 0);
    RNG ran2(id2, 1);

    // simple random check
    double r1 = 0.0, r2 = 0.0;
    for (int i = 0; i < 10000; i++)
    {
        r1 += ran1.ran();
        r2 += ran2.ran();
    }

    cout << endl;
    cout << "Simple randomness check, r = " << fixed << setw(12)
         << setprecision(6) << r1 / 10000.0 << " using 10000 samples (0.5)"
         << endl;
    cout << endl;

    EXPECT_SOFTEQ(0.5, r1/10000, 0.01);
    EXPECT_SOFTEQ(0.5, r2/10000, 0.01);

    RNG ran3(id3, 0);

    double r3 = 0.0;
    for (int i = 0; i < 10000; i++)
    {
        r3 += ran3.ran();
    }

    EXPECT_SOFTEQ(r1, r3, 1.0e-12);

    // now check that 1, 3 give equal streams
    double eps = 0.00001;
    for (int i = 0; i < 10; i++)
    {
        EXPECT_SOFTEQ(ran1.ran(), ran3.ran(), eps);
    }
}

//---------------------------------------------------------------------------//

TEST(RNG, Interface)
{
    int num  = 5;

    // make two sprng states
    int *idr = init_sprng(0, num, seed, 1);
    int *id1 = init_sprng(0, num, seed, 1);

    // now make some sprngs
    RNG ranr(idr, 0);
    RNG ran1(id1, 0);
    RNG ran2(ran1);
    RNG ran3(ran1);
    RNG empty;

    EXPECT_EQ(0, ran1.get_num());
    EXPECT_TRUE(ran1.assigned());
    EXPECT_TRUE(!empty.assigned());

    // get some reference numbers
    vector<double> ref(80);
    for (int i = 0; i < 80; i++)
    {
        ref[i] = ranr.ran();
    }

    // now check these sprngs that ALL access the same random number state
    double eps = .00001;
    for (int i = 0; i < 20; i++)
    {
        EXPECT_SOFTEQ(ref[i], ran1.ran(), eps);
    }
    for (int i = 20; i < 40; i++)
    {
        EXPECT_SOFTEQ(ref[i], ran2.ran(), eps);
    }
    for (int i = 40; i < 60; i++)
    {
        EXPECT_SOFTEQ(ref[i], ran3.ran(), eps);
    }

    // now check the ids
    EXPECT_EQ(id1, ran1.get_id());
    EXPECT_EQ(id1, ran2.get_id());
    EXPECT_EQ(id1, ran3.get_id());

    // now assignment
    ranr = ran2;
    for (int i = 60; i < 80; i++)
    {
        EXPECT_SOFTEQ(ref[i], ranr.ran(), eps);
    }
    EXPECT_EQ(id1, ranr.get_id());

    // check empty
    ran3 = empty;
    EXPECT_TRUE(!ran3.assigned());

    RNG ran4(ran3);
    EXPECT_TRUE(!ran4.assigned());
}

//---------------------------------------------------------------------------//

TEST(RNG, Pack)
{
    int num = 5;

    vector<double> ref(80);
    vector<char>   packed;

    {
        // make a sprng state
        int *id1 = init_sprng(0, num, seed, 1);
        int *idr = init_sprng(0, num, seed, 1);

        // now make some sprngs
        RNG ran1(id1, 0);
        RNG ranr(idr, 0);

        // get some reference numbers
        for (int i = 0; i < 80; i++)
            ref[i] = ranr.ran();

        // get 40 numbers from the non-ref
        for (int i = 0; i < 40; i++)
            ran1.ran();

        // pack up the sprng
        packed = ran1.pack();

        EXPECT_EQ(ran1.get_size(), packed.size());
    }

    // now check it
    {
        RNG uran(packed);

        // now check the first 40 Unpacked ran numbers
        double r   = 0;
        double rf  = 0;
        for (int i = 0; i < 40; i++)
        {
            r  = uran.ran();
            rf = ref[i+40];

            EXPECT_SOFTEQ(rf, r, 1.0e-12);
        }
    }
}

//---------------------------------------------------------------------------//
//                 end of tstRNG.cc
//---------------------------------------------------------------------------//
