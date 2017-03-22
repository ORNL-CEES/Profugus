//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CudaUtils/test/Device_View_Tester.hh
 * \author Tom Evans
 * \date   Fri Nov 18 11:55:09 2016
 * \brief  Device_View_Tester class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CudaUtils_test_Device_View_Tester_hh
#define CudaUtils_test_Device_View_Tester_hh

#include <string>
#include <vector>
#include <memory>

namespace host
{

//---------------------------------------------------------------------------//

class Monkey
{
  private:

    double d_height;
    int    d_id;

  public:

    Monkey(int    id,
           double height)
        : d_height(height)
        , d_id(id)
    {
    }

    int id() const { return d_id; }
    double height() const { return d_height; }
};

//---------------------------------------------------------------------------//

class Cage
{
  public:

    using Monkeys = std::vector<std::shared_ptr<Monkey>>;

  private:

    Monkeys d_monkeys;

  public:

    Cage(Monkeys monkeys)
        : d_monkeys(std::move(monkeys))
    {
    }

    const Monkeys& monkeys() const { return d_monkeys; }
};

//---------------------------------------------------------------------------//

class Zoo
{
  public:

    using Cages = std::vector<std::shared_ptr<Cage>>;

  private:

    Cages d_cages;

  public:
    Zoo(Cages cages)
        : d_cages(std::move(cages))
    {
    }

    const Cages& cages() const { return d_cages; }
};

}

//===========================================================================//
// TEST FUNCTIONS FOR DEVICE VIEW
//===========================================================================//

void test_one_level(const host::Zoo &zoo);
void test_two_level(const std::shared_ptr<host::Zoo> &zoo);

#endif // CudaUtils_test_Device_View_Tester_hh

//---------------------------------------------------------------------------//
// end of CudaUtils/test/Device_View_Tester.hh
//---------------------------------------------------------------------------//
