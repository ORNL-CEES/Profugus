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
#include <cmath>
#include <utility>
#include <iomanip>

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

namespace source
{

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

} // end namespace source

TEST(OMP, sp_copyin)
{
    using namespace source;

    omp_set_dynamic(0);
    omp_set_num_threads(4);

    // Build source on master thread
    Thread_Data::d_source    = std::make_shared<Source>();
    Thread_Data::d_source->N = 100;

    Transporter t;
    t.set_source();
}

//---------------------------------------------------------------------------//

namespace tally
{

enum Type {A, B, C};

class Tally
{
  private:
    int d_type;

    // Local tally
    double d_l;

    // Moments.
    double d_first;
    double d_second;

  public:
    Tally(Type type)
        : d_type(type)
        , d_l(0.0)
        , d_first(0.0)
        , d_second(0.0)
    {
        /* * */
    }

    void accumulate(Type t, double v)
    {
        if (t == d_type)
        {
            d_l += v;
        }
    }

    void end_history()
    {
        // Store moments
        d_first  += d_l;
        d_second += d_l * d_l;

        // Reset local
        d_l = 0.0;
    }

    void finalize(int N)
    {
        double inv_N = 1.0 / N;

        double l1 = d_first  * inv_N;
        double l2 = d_second * inv_N;

        double var = N / (static_cast<double>(N) - 1) * (l2 - l1*l1);

        d_first  = l1;
        d_second = std::sqrt(var / N);
    }

    double mean()  const { return d_first; }
    double error() const { return d_second; }
    int    type()  const { return d_type; }

    std::pair<double, double> get_moments() const
    {
        return std::make_pair(d_first, d_second);
    }

    void set_moments(std::pair<double, double> moments)
    {
        d_first  = moments.first;
        d_second = moments.second;
    }
};

class Tallier
{
  public:
    typedef std::shared_ptr<Tally> SP_Tally;
    typedef std::vector<SP_Tally>  Tallies;

  private:
    Tallies d_tallies;

  public:
    Tallier() { /* * */ }

    int num_tallies() const { return d_tallies.size(); }

    void add_tally(SP_Tally t)
    {
        d_tallies.push_back(t);
    }

    const Tallies& get_tallies() const
    {
        return d_tallies;
    }

    void accumulate(Type type, double v)
    {
        for (const auto &t : d_tallies)
        {
            t->accumulate(type, v);
        }
    }

    void end_history()
    {
        for (const auto &t : d_tallies)
        {
            t->end_history();
        }
    }

    void finalize(int N)
    {
        for (const auto &t : d_tallies)
        {
            t->finalize(N);
        }
    }
};

class Thread_Data
{
  public:
    static std::shared_ptr<Tallier> t;
#pragma omp threadprivate(t)
};

std::shared_ptr<Tallier> Thread_Data::t;

void build()
{
    omp_set_dynamic(0);
    omp_set_num_threads(4);

    // Build the talliers
#pragma omp parallel
    {
        Thread_Data::t = std::make_shared<Tallier>();

        // Make tallies and add them
        auto a = std::make_shared<Tally>(A);
        auto b = std::make_shared<Tally>(B);

        Thread_Data::t->add_tally(a);
        Thread_Data::t->add_tally(b);
    }
}

void thread_finalize(int N)
{
    REQUIRE(!profugus::in_thread_parallel_region());

    std::vector<double> l1(Thread_Data::t->num_tallies(), 0.0);
    std::vector<double> l2(Thread_Data::t->num_tallies(), 0.0);

#pragma omp parallel
    {
        const auto &tallies = Thread_Data::t->get_tallies();

        int tctr = 0;
        for (auto &t : tallies)
        {
            auto moments = t->get_moments();
#pragma omp atomic update
            l1[tctr] += moments.first;
#pragma omp atomic update
            l2[tctr] += moments.second;

            ++tctr;
        }
    }

    auto &tallies = Thread_Data::t->get_tallies();
    int tctr = 0;
    for (auto &t : tallies)
    {
        t->set_moments(std::make_pair(l1[tctr], l2[tctr]));
        ++tctr;
    }

    Thread_Data::t->finalize(N);
}

} // end namespace tally

TEST(OMP, threadprivate_tally)
{
    using namespace tally;

    build();

    // Do simulated sampling
    std::vector< std::vector< std::pair<Type, double> > > r = {
        {{A, 1.0}, {B, 2.0}, {B, 1.0}, {A, 3.0}, {C, 4.0}, {C, 5.0}},
        {{A, 2.0}, {B, 3.0}, {B, 5.0}, {A, 4.0}, {C, 9.0}, {C, 1.0}},
        {{A, 3.0}, {B, 4.0}, {B, 6.0}, {A, 2.0}, {C, 2.0}, {C, 9.0}},
        {{A, 4.0}, {B, 5.0}, {B, 7.0}, {A, 7.0}, {C, 1.0}, {C, 2.0}}};

    int N = 12;

#pragma omp parallel
    {
        int id = profugus::thread_id();

        double off = 0.0;

#pragma omp for
        for (int n = 0; n < N; ++n)
        {
            for (const auto &s : r[id])
            {
                Thread_Data::t->accumulate(s.first, s.second + off);
            }

            Thread_Data::t->end_history();

            off += 0.1;
        }
    }

    thread_finalize(N);

    const auto &tallies = Thread_Data::t->get_tallies();

    for (const auto &t : tallies)
    {
        if (t->type() == A)
        {
            EXPECT_FLOAT_EQ(6.7, t->mean());
            EXPECT_SOFTEQ(0.81333582, t->error(), 1.0e-6);
        }
        else if (t->type() == B)
        {
            EXPECT_FLOAT_EQ(8.45, t->mean());
            EXPECT_SOFTEQ(1.00968792, t->error(), 1.0e-6);
        }
        else
        {
            EXPECT_TRUE(0);
        }
    }
}

//---------------------------------------------------------------------------//
// end of MC/mc/test/tstOMP.cc
//---------------------------------------------------------------------------//
