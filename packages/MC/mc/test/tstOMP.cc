//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/test/tstOMP.cc
 * \author Thomas M. Evans
 * \date   Tue Aug 18 13:11:11 2015
 * \brief  OMP class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>
#include <vector>
#include <algorithm>

#include <omp.h>

#include "Utils/gtest/utils_gtest.hh"
#include "Utils/rng/RNG_Control.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class TestRNG : public testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::RNG_Control          RNG_Control_t;
    typedef RNG_Control_t::RNG_t           RNG_t;
    typedef std::shared_ptr<RNG_Control_t> SP_RNG_Control;
    typedef std::shared_ptr<RNG_t>         SP_RNG;

  protected:
    void SetUp()
    {
        rng_control = std::make_shared<RNG_Control_t>(324234);
        omp_set_dynamic(0);

        node  = profugus::node();
        nodes = profugus::nodes();
    }

    void build(int num_threads)
    {
        double begin = omp_get_wtime();

        omp_set_num_threads(num_threads);

        int total = nodes * omp_get_max_threads();
        cout << "Total RNGs = " << total << endl;

#pragma omp parallel
        {
#pragma omp critical
            {
                thread_id =
                    omp_get_thread_num() + omp_get_max_threads() * node;
                rng = rng_control->rng(thread_id);
                cout << "Building RNG for (" << omp_get_thread_num()
                     << "," << node << "); id = " << thread_id << endl;
            }
        }

        double end = omp_get_wtime();

        profugus::global_barrier();
    }

  protected:
    // >>> DATA

    SP_RNG_Control rng_control;

    int node, nodes;

    // >>> THREAD DATA
    static RNG_t rng;
#pragma omp threadprivate(rng)

    static int thread_id;
#pragma omp threadprivate(thread_id)
};

int TestRNG::thread_id = 0;
TestRNG::RNG_t TestRNG::rng;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(TestRNG, 1_thread)
{
    build(1);

#pragma omp parallel
    {
        int id = 0;
#pragma omp critical
        {
            auto refr = rng_control->rng(id);
            for (int n = 0; n < 100; ++n)
            {
                EXPECT_SOFTEQ(refr.ran(), rng.ran(), 1.0e-6);
            }
            ++id;
        }
    }
}

//---------------------------------------------------------------------------//

TEST_F(TestRNG, 4_thread)
{
    build(4);

#pragma omp parallel
    {
#pragma omp critical
        {
            cout << "Building reference rng for thread = " << thread_id
                 << endl;
            auto refr = rng_control->rng(thread_id);
            for (int n = 0; n < 1; ++n)
            {
                EXPECT_SOFTEQ(refr.ran(), rng.ran(), 1.0e-6);
            }
        }
    }
}

//---------------------------------------------------------------------------//

static int x[] = {0,1,2,3};
#pragma omp threadprivate(x)

TEST(OMP, threadprivate)
{
    omp_set_dynamic(0);
    omp_set_num_threads(4);

#pragma omp parallel
    {
        int id = omp_get_thread_num();
        x[id]  = 10;

#pragma omp critical
        {
            cout << "Values on " << id << " are "
                 << x[0] << "," << x[1] << "," << x[2] << "," << x[3]
                 << endl;
        }
    }
}

//---------------------------------------------------------------------------//

TEST(OMP, copyin)
{
    omp_set_dynamic(0);
    omp_set_num_threads(4);

    std::fill(std::begin(x), std::end(x), 0);

#pragma omp parallel copyin(x)
    {
        int id = omp_get_thread_num();
        x[id]  = 10;

#pragma omp critical
        {
            cout << "Values on " << id << " are "
                 << x[0] << "," << x[1] << "," << x[2] << "," << x[3]
                 << endl;
        }
    }
}

//---------------------------------------------------------------------------//

class Source
{
  public:
    double x, y, z;
    int N;
};

class Thread_Data
{
  public:

    static std::shared_ptr<Source> d_source;
#pragma omp threadprivate(d_source)
};

std::shared_ptr<Source> Thread_Data::d_source;

class Transporter
{
  private:
    std::shared_ptr<Source> d_source;

  public:
    Transporter()
    {
        /* * */
    }

    void set_source()
    {
#pragma omp parallel copyin(Thread_Data::d_source)
        {
            int id   = omp_get_thread_num();

#pragma omp critical
            {
                cout << "Thread " << id << " has " << Thread_Data::d_source->N
                     << " particles" << endl;
            }

            d_source = Thread_Data::d_source;
            int nt   = omp_get_num_threads();

            d_source->N /= nt;

#pragma omp critical
            {
                cout << "Thread " << id << " has " << d_source->N
                     << " particles" << endl;
            }
        }
    }
};

TEST(OMP, sp_copyin)
{
    omp_set_dynamic(0);
    omp_set_num_threads(4);

    // Build source on master thread
    Thread_Data::d_source    = std::make_shared<Source>();
    Thread_Data::d_source->N = 100;

    Transporter t;
    t.set_source();
}

//---------------------------------------------------------------------------//
// end of MC/mc/test/tstOMP.cc
//---------------------------------------------------------------------------//
