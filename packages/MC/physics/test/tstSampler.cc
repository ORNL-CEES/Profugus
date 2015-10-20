//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstSampler.cc
 * \author Thomas M. Evans
 * \date   Friday May 2 10:33:33 2014
 * \brief  Sampler functions test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Sampler.hh"

#include <vector>

#include "gtest/utils_gtest.hh"
#include "utils/Vector_Functions.hh"
#include "rng/RNG.hh"


using namespace std;

using profugus::sampler::sample_discrete_CDF;
using profugus::sampler::sample_dcdf;
using profugus::sampler::sample_small_dcdf;
using profugus::sampler::sample_smallrev_dcdf;
using profugus::sampler::sample_watt;
using profugus::sampler::sample_linear;
using profugus::RNG;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(DiscreteCDF, MultiBin)
{
    // make a CDF
    const double cdf[] = { 0.21, 0.32, 0.64, 0.99, 1.00};

    // boundary samples
    EXPECT_EQ(0, sample_discrete_CDF(5, cdf, 0.21));
    EXPECT_EQ(1, sample_discrete_CDF(5, cdf, 0.32));
    EXPECT_EQ(2, sample_discrete_CDF(5, cdf, 0.64));
    EXPECT_EQ(3, sample_discrete_CDF(5, cdf, 0.99));
    EXPECT_EQ(4, sample_discrete_CDF(5, cdf, 1.00));

    // test some samples
    EXPECT_EQ(0, sample_discrete_CDF(5, cdf, 0.0000000));
    EXPECT_EQ(0, sample_discrete_CDF(5, cdf, 0.2099999));
    EXPECT_EQ(1, sample_discrete_CDF(5, cdf, 0.2100001));
    EXPECT_EQ(2, sample_discrete_CDF(5, cdf, 0.3200001));
    EXPECT_EQ(2, sample_discrete_CDF(5, cdf, 0.6399999));
    EXPECT_EQ(3, sample_discrete_CDF(5, cdf, 0.6400001));
    EXPECT_EQ(3, sample_discrete_CDF(5, cdf, 0.9899999));
    EXPECT_EQ(4, sample_discrete_CDF(5, cdf, 0.9900001));
}

TEST(DiscreteCDF, OneBin)
{
    // 1-bin cdf
    const double cdf[] = {1.00};
    EXPECT_EQ(0, sample_discrete_CDF(1, cdf, 0.0));
    EXPECT_EQ(0, sample_discrete_CDF(1, cdf, 0.2));
    EXPECT_EQ(0, sample_discrete_CDF(1, cdf, 0.3));
    EXPECT_EQ(0, sample_discrete_CDF(1, cdf, 0.9));
    EXPECT_EQ(0, sample_discrete_CDF(1, cdf, 1.0));
}

TEST(DiscreteCDF, ZeroProbabilities)
{
    // make a CDF with some zero probabilities
    const double cdf[] = {0.21, 0.32, 0.32, 0.32, 1.00};

    EXPECT_EQ(0, sample_discrete_CDF(5, &cdf[0], 0.209));
    EXPECT_EQ(1, sample_discrete_CDF(5, &cdf[0], 0.211));
    EXPECT_EQ(1, sample_discrete_CDF(5, &cdf[0], 0.319));
    EXPECT_EQ(4, sample_discrete_CDF(5, &cdf[0], 0.321));
}

TEST(SampleDCDF, MultiBin)
{
    // make a CDF
    const double cdf[] = { 0.21, 0.32, 0.64, 0.99, 1.00};

    EXPECT_EQ(0, sample_dcdf(0.21        , cdf, std::end(cdf)) - cdf);
    EXPECT_EQ(1, sample_dcdf(0.32        , cdf, std::end(cdf)) - cdf);
    EXPECT_EQ(4, sample_dcdf(0.9999999999, cdf, std::end(cdf)) - cdf);

    EXPECT_EQ(0, sample_dcdf(0.0      , cdf, std::end(cdf)) - cdf);
    EXPECT_EQ(3, sample_dcdf(0.6400001, cdf, std::end(cdf)) - cdf);
    EXPECT_EQ(4, sample_dcdf(0.999    , cdf, std::end(cdf)) - cdf);
}

TEST(SampleSmallDCDF, MultiBin)
{
    // make a CDF
    const double cdf[] = {0.21, 0.32, 0.64, 0.99, 1.00};

    EXPECT_EQ(0, sample_small_dcdf(0.21        , cdf, std::end(cdf)) - cdf);
    EXPECT_EQ(1, sample_small_dcdf(0.32        , cdf, std::end(cdf)) - cdf);
    EXPECT_EQ(4, sample_small_dcdf(0.9999999999, cdf, std::end(cdf)) - cdf);

    EXPECT_EQ(0, sample_small_dcdf(0.0      , cdf, std::end(cdf)) - cdf);
    EXPECT_EQ(3, sample_small_dcdf(0.6400001, cdf, std::end(cdf)) - cdf);
    EXPECT_EQ(4, sample_small_dcdf(0.999    , cdf, std::end(cdf)) - cdf);
}

TEST(SampleSmallRevDCDF, MultiBin)
{
    // make a CDF
    const double cdf[] = {0.21, 0.32, 0.64, 0.99, 1.00};

    EXPECT_EQ(0, sample_smallrev_dcdf(0.21        , cdf, std::end(cdf)) - cdf);
    EXPECT_EQ(1, sample_smallrev_dcdf(0.32        , cdf, std::end(cdf)) - cdf);
    EXPECT_EQ(4, sample_smallrev_dcdf(0.9999999999, cdf, std::end(cdf)) - cdf);

    EXPECT_EQ(0, sample_smallrev_dcdf(0.0      , cdf, std::end(cdf)) - cdf);
    EXPECT_EQ(3, sample_smallrev_dcdf(0.6400001, cdf, std::end(cdf)) - cdf);
    EXPECT_EQ(4, sample_smallrev_dcdf(0.999    , cdf, std::end(cdf)) - cdf);
}

//---------------------------------------------------------------------------//

TEST(Watt, Samples)
{
    const int N = 25000;  // number of samples
    int tally[4] = {0};   // energy bin scores
    double sum = 0.0;     // sum of sampled energies

    // make a SPRNG object
    RNG rng(0);

    // take samples, sum, tally
    for (int i = 0; i < N; i++)
    {
        double E = sample_watt(rng);
        EXPECT_GT(E, 0);

        sum += E;

        // bin tallies into 4 equi-probable energy ranges
        if (E < 0.8298e6)
            ++tally[0];
        else if (E < 1.596e6)
            ++tally[1];
        else if (E < 2.722e6)
            ++tally[2];
        else
            ++tally[3];
    }

    // test mean energy value (~ 1.98 MeV)
    EXPECT_SOFTEQ(sum / N, 1.98e6, 0.01);

    // test energy bin tallies
    for (int i = 0; i < 4; i++)
    {
        EXPECT_SOFTEQ(tally[i] / static_cast<double>(N), 0.25, .1)
            << "failure in bin " << i;
    }
}

//---------------------------------------------------------------------------//

TEST(SampleLinear, Samples)
{
    const int N = 25000;  // number of samples
    double sum = 0.0;     // sum of sampled energies
    double mean;

    // make a RNG object
    RNG rng(0);

    // Loop over samples
    for(int i = 0; i < N; ++i)
    {
        sum += sample_linear(rng.ran());
    }

    // Calc mean (exact = 2/3)
    mean = sum / static_cast<double>(N);
    double eps = 0.05;

    EXPECT_TRUE(mean >= (2.0/3.0)-eps && mean <= (2.0/3.0)+eps);
}

//---------------------------------------------------------------------------//

TEST(SampleArbLinear, Samples)
{
    const int N = 100000;  // number of samples
    double sum = 0.0;     // sum of sampled energies
    double mean;

    // make a RNG object
    RNG rng(0);

    // Loop over samples
    for(int i = 0; i < N; ++i)
    {
        sum += sample_linear<double>(rng, 1.0, 2.0);
    }

    // Calc mean (exact = 5/9)
    mean = sum / static_cast<double>(N);
    double eps = std::sqrt(N);

    EXPECT_SOFTEQ(5.0/9.0, mean, eps);
}

//---------------------------------------------------------------------------//
//                        end of tstSampler.cc
//---------------------------------------------------------------------------//
