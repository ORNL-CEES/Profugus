//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstStream.cc
 * \author Seth R Johnson
 * \date   Wed Oct 02 15:11:19 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../cuda_utils/Stream.hh"

#include <utility>
#include "gtest/utils_gtest.hh"

#include "../cuda_utils/Device_Vector.hh"
#include "../cuda_utils/Host_Vector.hh"
#include "../cuda_utils/Hardware.hh"
#include "../cuda_utils/Event.hh"

#include "Stream_Test_Kernel.cuh"
#include "Stream_Test_Kernel_Data.hh"

using cuda::Stream;

#if DBC == 0
#define NUM_SEGMENTS 128u
#else
#define NUM_SEGMENTS 16u
#endif

//---------------------------------------------------------------------------//
TEST(StreamUsageTest, copy_constructor)
{
    typedef Stream<cuda::arch::Host> Stream_t;

    Stream_t s;
    EXPECT_EQ(1, s.use_count());
    {
        // Copy constructor
        Stream_t s2(s);
        EXPECT_EQ(2, s.use_count());
        EXPECT_EQ(2, s2.use_count());
    }
    EXPECT_EQ(1, s.use_count());
}

//---------------------------------------------------------------------------//
TEST(StreamUsageTest, swap)
{
    typedef Stream<cuda::arch::Host> Stream_t;

    Stream_t sa;
    Stream_t sa2(sa);
    Stream_t sb;

    EXPECT_EQ(2, sa.use_count());
    EXPECT_EQ(2, sa2.use_count());
    EXPECT_EQ(1, sb.use_count());

    swap(sa, sb);
    EXPECT_EQ(1, sa.use_count());
    EXPECT_EQ(2, sa2.use_count());
    EXPECT_EQ(2, sb.use_count());
}

//---------------------------------------------------------------------------//
TEST(StreamUsageTest, assign)
{
    typedef Stream<cuda::arch::Host> Stream_t;

    Stream_t sa;
    Stream_t sa2(sa);
    Stream_t sb;

    EXPECT_EQ(2, sa.use_count());
    EXPECT_EQ(2, sa2.use_count());
    EXPECT_EQ(1, sb.use_count());

    sb = sa;
    EXPECT_EQ(3, sa.use_count());
    EXPECT_EQ(3, sa2.use_count());
    EXPECT_EQ(3, sb.use_count());
}

//---------------------------------------------------------------------------//
template<typename Arch_T>
struct Completion
{
    enum { HAS_ASYNC = 0 };
};

template<>
struct Completion<cuda::arch::Device>
{
    enum { HAS_ASYNC = 1 };
};

//---------------------------------------------------------------------------//
template<typename Arch_T, typename T>
struct Kernel_Traits
{
    typedef Arch_T Arch_t;
    typedef T      Float_t;
};

//---------------------------------------------------------------------------//
// POLYGLOT TESTS
//---------------------------------------------------------------------------//
template<typename Kernel_Traits>
class StreamTest : public ::testing::Test
{
  protected:
    typedef Kernel_Traits                     Kernel_Traits_t;
    typedef typename Kernel_Traits_t::Arch_t  Arch_t;
    typedef typename Kernel_Traits_t::Float_t Float_t;

    typedef cuda::Hardware<Arch_t>               Hardware_t;
    typedef cuda::Event<Arch_t>                  Event_t;
    typedef cuda::Stream<Arch_t>                 Stream_t;
    typedef cuda::Device_Vector<Arch_t, Float_t> Device_Vector_t;
    typedef cuda::Host_Vector<Float_t>           Host_Vector_t;

    typedef cuda::Stream_Test_Kernel_Data<Arch_t, Float_t> Kernel_Data_t;

    void SetUp()
    {
        // Initialize device
        if (!Hardware_t::have_acquired())
        {
            std::cout << "Acquiring device..." << std::flush;
            Hardware_t::acquire();
            std::cout << "done." << std::endl;
        }
        INSIST(Hardware_t::have_acquired(), "Device could not be acquired.");

        d_num_blocks   = 8 * Hardware_t::num_multiprocessors();
        d_num_threads  = 2 * Hardware_t::num_cores_per_mp();
        d_num_segments = NUM_SEGMENTS;
    }

    void build_local_tau(Host_Vector_t& tau, Float_t scale) const
    {
        REQUIRE(tau.size() == num_rays() * d_num_segments);
        REQUIRE(scale > 0. && scale < 1.);

        typename Host_Vector_t::iterator it = tau.begin();

        for (unsigned int seg = 0; seg < d_num_segments; ++seg)
        {
            for (unsigned int ray = 0, ray_end = num_rays();
                    ray != ray_end; ++ray, ++it)
            {
                *it = static_cast<Float_t>(1) / ((seg + ray) % 4 + 2);
                *it *= scale;
            }
        }
        CHECK(it == tau.end());
    }

    void build_local_src(Host_Vector_t& src) const
    {
        REQUIRE(src.size() == num_rays() * d_num_segments);

        typename Host_Vector_t::iterator it = src.begin();

        for (unsigned int seg = 0; seg < d_num_segments; ++seg)
        {
            for (unsigned int ray = 0, ray_end = num_rays();
                    ray != ray_end; ++ray, ++it)
            {
                *it = (seg + ray) % 3;
            }
        }
        CHECK(it == src.end());
    }

    unsigned int num_rays() const
    {
        return std::max(d_num_blocks * d_num_threads, 32u);
    }

  protected:
    // Number of blocks to use
    unsigned int d_num_blocks;
    // Number of threads to use
    unsigned int d_num_threads;
    // Number of segments per thread
    std::size_t  d_num_segments;
};


typedef Kernel_Traits<cuda::arch::Host, float>  KT_HF;
typedef Kernel_Traits<cuda::arch::Host, double> KT_HD;
#ifdef USE_CUDA
typedef Kernel_Traits<cuda::arch::Device, float>  KT_DF;
typedef Kernel_Traits<cuda::arch::Device, double> KT_DD;
// instantiate both host and device code
typedef ::testing::Types<KT_HF, KT_HD, KT_DF, KT_DD> ArchDataTypes;
//typedef ::testing::Types<KT_DD> ArchDataTypes;
#else
// instantiate host-only code
typedef ::testing::Types<KT_HF, KT_HD> ArchDataTypes;
#endif

TYPED_TEST_CASE(StreamTest, ArchDataTypes);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TYPED_TEST(StreamTest, performance)
{
    typedef typename TestFixture::Arch_t        Arch_t;
    typedef typename TestFixture::Event_t       Event_t;
    typedef typename TestFixture::Host_Vector_t Host_Vector_t;
    typedef typename TestFixture::Kernel_Data_t Kernel_Data_t;
    typedef typename TestFixture::Stream_t      Stream_t;

    // Events for keeping track of when the kernels were started and stopped
    std::cout << "Creating events..." << std::endl;
    Event_t copy_stop, kernel_start, kernel_stop;

    // Create front and back streams
    std::cout << "Creating streams..." << std::endl;
    Stream_t active_stream;
    Stream_t build_stream;

    EXPECT_EQ(1, active_stream.use_count());
    EXPECT_EQ(1, build_stream.use_count());

    // Initialize kernel data and launch arguments
    // Each kernel data instance (launch args) gets its own stream.
    Kernel_Data_t data(this->num_rays(), this->d_num_segments);
    data.launch_args.block_size = this->d_num_threads;
    data.launch_args.grid_size  = this->d_num_blocks;
    data.launch_args.stream     = active_stream;

    // Whether we can expect asynchronous behavior
    const bool complete = (Completion<Arch_t>::HAS_ASYNC ? false : true);

    // Host vector for building tau
    Host_Vector_t host_tau(data.num_rays * data.num_segments);
    Host_Vector_t host_src(data.num_rays * data.num_segments);
    Host_Vector_t host_input(data.num_rays);
    Host_Vector_t host_output(data.num_rays);

    std::cout << "Running on " << data.launch_args.grid_size << " blocks"
        << " x " << data.launch_args.block_size << " threads..." << std::flush;

    // First initialization of tau
    std::cout << "T" << std::flush;
    this->build_local_tau(host_tau, 1 / (2.));
    // Asynchronously copy tau
    std::cout << ">" << std::flush;
    data.tau_build.assign_async(host_tau, build_stream);

    for (unsigned int i = 0; i < 10; ++i)
    {
        EXPECT_TRUE(active_stream.is_complete());
        // Always rebuild source
        std::cout << "S" << std::flush;
        this->build_local_src(host_src);
        // Asynchronously copy src
        std::cout << ">" << std::flush;
        data.source.assign_async(host_src, build_stream);
        // Event lets us know when the copy is complete
        copy_stop.record(build_stream);
        EXPECT_EQ(complete, copy_stop.is_complete());

        // Always rebuild input
        std::cout << "I" << std::flush;
        std::fill(host_input.begin(), host_input.end(), 1.);
        // Asynchronously copy input
        std::cout << ">" << std::flush;
        data.input.assign_async(host_input, build_stream);

        // Wait until building/copy of the prior tau is completed
        build_stream.synchronize();
        EXPECT_TRUE(build_stream.is_complete());

        // Swap in the now-built tau (pointer swap, cheap)
        swap(data.tau_build, data.tau);

        // Call the kernel
        std::cout << "K" << std::flush;
        kernel_start.record(active_stream);
        cuda::stream_test(data);
        kernel_stop.record(active_stream);
        EXPECT_EQ(complete, active_stream.is_complete());

        // precompute next tau while kernel executes
        std::cout << "T" << std::flush;
        this->build_local_tau(host_tau, 1 / (i + 3.));
        // Asynchronously copy tau
        std::cout << ">" << std::flush;
        data.tau_build.assign_async(host_tau, build_stream);
        EXPECT_EQ(complete, build_stream.is_complete());

        // Asynchronously copy solution to host after kernel completes
        std::cout << "<" << std::flush;
        device_to_host_async(data.output, host_output, active_stream);
        // Copy from kernel should still be going
        EXPECT_EQ(complete, active_stream.is_complete());

        // Wait until the kernel and the "stop" event are completed
        kernel_stop.synchronize();
        // But now the asynchronous call to device_to_host will probably be
        // going
        EXPECT_EQ(complete, active_stream.is_complete());
        std::cout << "!" << std::flush;
        active_stream.synchronize();
        // The synchronize to the event should ensure that the stream finishes
        EXPECT_TRUE(active_stream.is_complete());
    }
    cout << "\n" << "Done." << endl;
}
//---------------------------------------------------------------------------//
//                 end of tstStream.cc
//---------------------------------------------------------------------------//
