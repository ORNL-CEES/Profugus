//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/test/tstEvent.cc
 * \author Seth R Johnson
 * \date   Thu Jul 11 09:47:38 2013
 * \brief
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Event.hh"

#include "Utils/gtest/utils_gtest.hh"

#include "../Device_Vector.hh"
#include "../Hardware.hh"
#include "../CudaDBC.hh"
#include "../Vector_Traits.hh"

#include "Profiler_Kernel.cuh"

#if DBC == 0
#define MIN_NUM_THREADS 512u
#else
#define MIN_NUM_THREADS 64u
#endif

//---------------------------------------------------------------------------//
template <typename Arch_Switch, typename T>
struct Kernel_Traits
{
    typedef Arch_Switch Arch_t;
    typedef T           Data_t;
};

//---------------------------------------------------------------------------//
// POLYGLOT TESTS
//---------------------------------------------------------------------------//
template<typename Kernel_Traits>
class EventTest : public ::testing::Test
{
  protected:
    typedef Kernel_Traits                    Kernel_Traits_t;
    typedef typename Kernel_Traits_t::Arch_t Arch_t;
    typedef typename Kernel_Traits_t::Data_t Data_t;

    typedef cuda::Hardware<Arch_t>              Hardware_t;
    typedef cuda::Event<Arch_t>                 Event_t;
    typedef cuda::Device_Vector<Arch_t, Data_t> Device_Vector_t;

    void SetUp()
    {
        // Initialize device
        if (!Hardware_t::have_acquired())
        {
            std::cout << "Acquiring device..." << std::flush;
            Hardware_t::acquire();
            std::cout << "done." << std::endl;
        }
        Insist(Hardware_t::have_acquired(), "Device could not be acquired.");

        d_num_blocks  = 8 * Hardware_t::num_multiprocessors();
        d_num_threads = 2 * Hardware_t::num_cores_per_mp();
        d_size        = std::max(d_num_blocks * d_num_threads, MIN_NUM_THREADS);
    }

  protected:
    // Number of blocks to use
    unsigned int d_num_blocks;
    // Number of threads to use
    unsigned int d_num_threads;
    // Number of elements to use
    std::size_t  d_size;
};


typedef Kernel_Traits<cuda::arch::Host, int>    KT_HI;
typedef Kernel_Traits<cuda::arch::Host, float>  KT_HF;
typedef Kernel_Traits<cuda::arch::Host, double> KT_HD;
#ifdef USE_CUDA
typedef Kernel_Traits<cuda::arch::Device, int>    KT_DI;
typedef Kernel_Traits<cuda::arch::Device, float>  KT_DF;
typedef Kernel_Traits<cuda::arch::Device, double> KT_DD;
// instantiate both host and device code
typedef ::testing::Types<KT_HI, KT_HF, KT_HD, KT_DI, KT_DF, KT_DD> ArchDataTypes;
#else
// instantiate host-only code
typedef ::testing::Types<KT_HI, KT_HF, KT_HD> ArchDataTypes;
#endif

TYPED_TEST_CASE(EventTest, ArchDataTypes);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TYPED_TEST(EventTest, performance)
{
    typedef typename TestFixture::Event_t         Event_t;
    typedef typename TestFixture::Device_Vector_t Device_Vector_t;

    const unsigned int num_blocks = this->d_num_blocks;
    const unsigned int num_threads = this->d_num_threads;

    // Accumulated run time *in milliseconds*
    double accum_time = 0.;
    unsigned long int operations = 0;

    // Events for keeping track of when the kernels were started and stopped
    std::cout << "Creating events..." << std::endl;
    Event_t event_start, event_stop;

    Device_Vector_t data(this->d_size);

    // Call the kernel once to "prime" CUDA (load dlls, generate binary
    // instructions for the current GPU choice)

    std::cout << "Doing first kernel executions..." << std::endl;
    cuda::operation_test(num_blocks, num_threads, data);

    std::cout << "Running on " << num_blocks << " blocks"
        << " x " << num_threads << " threads..." << std::flush;

    for (unsigned int i = 0; i != 20; ++i)
    {
        // Start the event before the kernel call
        event_start.record();

        // Call the CUDA code which wraps the kernel
        operations += cuda::operation_test(num_blocks, num_threads, data);

        // Add an event that fires after the kernel is finished
        event_stop.record();

        // Wait until the kernel and the "stop" event are completed
        event_stop.synchronize();
        // Calculate the time delta between start and stop
        float this_time = event_stop.elapsed_time_since(event_start);

        // Make sure the time is positive if we've run sufficient samples
        if (num_blocks * num_threads > 256)
            EXPECT_GT(this_time, 0.f);
        // Accumulate milliseconds
        accum_time += this_time;

        cout << "." << std::flush;
    }
    cout << std::endl;

    cout << "Performance (Gflops): "
         << (operations / (accum_time / 1000)) / (1000ul * 1000 * 1000)
         << endl;
}

//---------------------------------------------------------------------------//
//                        end of tstEvent.cc
//---------------------------------------------------------------------------//
