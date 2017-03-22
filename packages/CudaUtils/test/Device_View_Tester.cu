//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   CudaUtils/test/Device_View_Tester.cu
 * \author Tom Evans
 * \date   Fri Nov 18 11:55:09 2016
 * \brief  Device_View_Tester kernel definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <memory>

#include "Device_View_Tester.hh"
#include "gtest/Gtest_Functions.hh"

#include "cuda_utils/Shared_Device_Ptr.hh"
#include "cuda_utils/Device_View.hh"

//---------------------------------------------------------------------------//
// DEVICE WRAPPERS
//---------------------------------------------------------------------------//

namespace device
{

class Monkey
{
  private:

    double d_height;
    int    d_id;

  public:

    Monkey(const std::shared_ptr<host::Monkey> &monkey)
        : d_height(monkey->height())
        , d_id(monkey->id())
    {
    }

    __device__
    int id() const { return d_id; }

    __device__
    double height() const { return d_height; }
};

//---------------------------------------------------------------------------//

class Cage
{
  public:

  private:

    Monkey *d_monkeys;
    int     d_size;

  public:

    Cage(Monkey *monkeys,
         int     size)
        : d_monkeys(monkeys)
        , d_size(size)
    {
    }

    __device__
    int size() const { return d_size; }

    __device__
    Monkey* monkeys() const { return d_monkeys; }
};

class Cage_Host_Manager : public cuda::Device_Memory_Manager<Cage>
{
  private:

    // Managed device memory
    thrust::device_vector<Monkey> d_monkeys;

    // Host object
    std::shared_ptr<host::Cage> d_host_cage;

  public:

    Cage_Host_Manager(const std::shared_ptr<host::Cage> &cage)
        : d_host_cage(cage)
    {
        std::vector<Monkey> monkeys;
        for (const auto &monkey : d_host_cage->monkeys())
        {
            monkeys.push_back(monkey);
        }
        d_monkeys = thrust::device_vector<Monkey>(
            monkeys.begin(), monkeys.end());
    }

    Cage device_instance()
    {
        Cage cage(d_monkeys.data().get(), d_monkeys.size());
        return cage;
    }
};

//---------------------------------------------------------------------------//

class Zoo
{
  private:

    Cage *d_cages;
    int   d_size;

  public:

    Zoo(Cage *cages,
        int   size)
        : d_cages(cages)
        , d_size(size)
    {
    }

    __device__
    int size() const { return d_size; }

    __device__
    Cage* cages() const { return d_cages; }
};

class Zoo_Host_Manager : public cuda::Device_Memory_Manager<Zoo>
{
  private:

    // Managed device memory
    std::shared_ptr<cuda::Device_View<Cage>> d_cages;

    // Host object
    std::shared_ptr<host::Zoo> d_host_zoo;

  public:

    Zoo_Host_Manager(const std::shared_ptr<host::Zoo> &zoo)
        : d_host_zoo(zoo)
    {
        using View     = cuda::Device_View<device::Cage>;
        using Managers = View::Managers;

        // Make a vector of cage managers
        Managers managers;
        for (const auto &cage : d_host_zoo->cages())
        {
            managers.push_back(
                std::make_shared<device::Cage_Host_Manager>(cage));
        }

        d_cages = std::make_shared<View>(managers);
    }

    Zoo device_instance()
    {
        Zoo zoo(d_cages->get_device_ptr(), d_cages->size());
        return zoo;
    }
};

}

//---------------------------------------------------------------------------//
// KERNELS
//---------------------------------------------------------------------------//

__global__
void one_level_kernel(const device::Cage *cages,
                      int                 size,
                      int                *ids)
{
    int ctr = 0;
    for (int n = 0; n < size; ++n)
    {
        auto monkeys = cages[n].monkeys();
        for (int m = 0; m < cages[n].size(); ++m)
        {
            ids[ctr] = monkeys[m].id();
            ++ctr;
        }
    }
}

//---------------------------------------------------------------------------//

__global__
void two_level_kernel(const device::Zoo *zoo,
                      int               *ids)
{
    auto cages = zoo->cages();
    int size   = zoo->size();

    int ctr = 0;
    for (int n = 0; n < size; ++n)
    {
        auto monkeys = cages[n].monkeys();
        for (int m = 0; m < cages[n].size(); ++m)
        {
            ids[ctr] = monkeys[m].id();
            ++ctr;
        }
    }
}

//---------------------------------------------------------------------------//
// CUDA TEST FUNCTIONS
//---------------------------------------------------------------------------//

void test_one_level(const host::Zoo &zoo)
{
    using View     = cuda::Device_View<device::Cage>;
    using Managers = View::Managers;

    thrust::device_vector<int> ids(5, -1);

    // Make a vector of cage managers
    Managers managers;
    for (const auto &cage : zoo.cages())
    {
        managers.push_back(std::make_shared<device::Cage_Host_Manager>(cage));
    }

    View view(managers);

    one_level_kernel<<<1,1>>>(view.get_device_ptr(),
                              view.size(),
                              ids.data().get());

    thrust::host_vector<int> result(ids.begin(), ids.end());
    EXPECT_VEC_EQ(result, std::vector<int>({1,2,3,4,5}));

}

//---------------------------------------------------------------------------//

void test_two_level(const std::shared_ptr<host::Zoo> &zoo)
{
    device::Zoo_Host_Manager zoo_manager(zoo);
    auto zoo_ptr = cuda_utils::shared_device_ptr<device::Zoo>(
        zoo_manager.device_instance());

    thrust::device_vector<int> ids(5, -1);

    two_level_kernel<<<1,1>>>(zoo_ptr.get_device_ptr(),
                              ids.data().get());

    thrust::host_vector<int> result(ids.begin(), ids.end());
    EXPECT_VEC_EQ(result, std::vector<int>({1,2,3,4,5}));
}

//---------------------------------------------------------------------------//
// end of CudaUtils/test/Device_View_Tester.cu
//---------------------------------------------------------------------------//
