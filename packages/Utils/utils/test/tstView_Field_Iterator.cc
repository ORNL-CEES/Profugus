//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/test/tstView_Field_Iterator.cc
 * \author Gregory G. Davidson
 * \date   Wed Oct 08 14:22:47 2014
 * \brief  Tests the View_Field_Iterator
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../View_Field_Iterator.hh"

#include "Utils/gtest/profugus_gtest.hh"

#include <vector>

using profugus::VF_Iterator;
using profugus::const_VF_Iterator;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class View_Field_IteratorTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef std::vector<int>       Vec_Int;
    typedef VF_Iterator<int>       VFI;
    typedef const_VF_Iterator<int> cVFI;

  protected:
    void SetUp()
    {
        int vals[] = {1,  2,  3,  4,
                      5,  6,  7,  8,
                      9,  10, 11, 12,
                      13, 14, 15, 16,
                      17, 18, 19, 20};
        vi = Vec_Int(std::begin(vals), std::end(vals));
    }

  protected:
    Vec_Int vi;
};

//---------------------------------------------------------------------------//
// VF_Iterator TESTS
//---------------------------------------------------------------------------//

TEST_F(View_Field_IteratorTest, stride_1_test)
{
    // Default constructor
    VFI def_vfi;
    EXPECT_EQ(1, def_vfi.stride());
    EXPECT_FALSE(def_vfi.is_valid());

    // Constructor
    VFI vfi(vi.data(), 1);
    EXPECT_EQ(1, vfi.stride());
    EXPECT_EQ(1, *(vfi.get_pointer()));
    EXPECT_TRUE(vfi.is_valid());

    // Copy constructor
    VFI vfi2(vfi);
    EXPECT_EQ(1, vfi2.stride());
    EXPECT_EQ(1, *(vfi2.get_pointer()));

    // Copy operator
    VFI vfi3;
    vfi3 = vfi2;
    EXPECT_EQ(1, vfi3.stride());
    EXPECT_EQ(1, *(vfi3.get_pointer()));

    // Dereference
    EXPECT_EQ(1, *vfi);
    // Prefix increment
    ++vfi;
    EXPECT_EQ(2, *vfi);
    // Prefix decrement
    --vfi;
    EXPECT_EQ(1, *vfi);

    // Postfix increment
    VFI vfi4 = vfi++;
    EXPECT_EQ(1, *vfi4);
    EXPECT_EQ(2, *vfi);
    // Postfix decrement
    VFI vfi5 = vfi--;
    EXPECT_EQ(2, *vfi5);
    EXPECT_EQ(1, *vfi);

    // Addition operator
    VFI vfi6 = vfi+2;
    EXPECT_EQ(3, *vfi6);
    // Subtraction operator
    VFI vfi7 = vfi6-2;
    EXPECT_EQ(1, *vfi7);
    EXPECT_EQ(2, vfi6 - vfi7);

    // Compound addition
    vfi += 2;
    EXPECT_EQ(3, *vfi);
    // Compound subtraction
    vfi -= 2;
    EXPECT_EQ(1, *vfi);

    // Comparators
    VFI vfi8 = vfi + 1;
    EXPECT_TRUE(vfi8 > vfi);
    EXPECT_TRUE(vfi8 >= vfi);
    EXPECT_TRUE(vfi == vfi);
    EXPECT_TRUE(vfi >= vfi);
    EXPECT_TRUE(vfi < vfi8);
    EXPECT_TRUE(vfi <= vfi8);

    // Offset
    EXPECT_EQ(3, vfi[2]);

    // Change the value
    *vfi = 100;
    EXPECT_EQ(100, *vfi);
}

//---------------------------------------------------------------------------//

TEST_F(View_Field_IteratorTest, stride_4_test)
{
    // Default constructor
    VFI def_vfi;
    EXPECT_EQ(1, def_vfi.stride());

    // Constructor
    VFI vfi(vi.data(), 4);
    EXPECT_EQ(4, vfi.stride());
    EXPECT_EQ(1, *(vfi.get_pointer()));

    // Copy constructor
    VFI vfi2(vfi);
    EXPECT_EQ(4, vfi2.stride());
    EXPECT_EQ(1, *(vfi2.get_pointer()));

    // Copy operator
    VFI vfi3;
    vfi3 = vfi2;
    EXPECT_EQ(4, vfi3.stride());
    EXPECT_EQ(1, *(vfi3.get_pointer()));

    // Dereference
    EXPECT_EQ(1, *vfi);
    // Prefix increment
    ++vfi;
    EXPECT_EQ(5, *vfi);
    // Prefix decrement
    --vfi;
    EXPECT_EQ(1, *vfi);

    // Postfix increment
    VFI vfi4 = vfi++;
    EXPECT_EQ(1, *vfi4);
    EXPECT_EQ(5, *vfi);
    // Postfix decrement
    VFI vfi5 = vfi--;
    EXPECT_EQ(5, *vfi5);
    EXPECT_EQ(1, *vfi);

    // Addition operator
    VFI vfi6 = vfi+2;
    EXPECT_EQ(9, *vfi6);
    // Subtraction operator
    VFI vfi7 = vfi6-2;
    EXPECT_EQ(1, *vfi7);
    EXPECT_EQ(2, vfi6 - vfi7);

    // Compound addition
    vfi += 2;
    EXPECT_EQ(9, *vfi);
    // Compound subtraction
    vfi -= 2;
    EXPECT_EQ(1, *vfi);

    // Comparators
    VFI vfi8 = vfi + 1;
    EXPECT_TRUE(vfi8 > vfi);
    EXPECT_TRUE(vfi8 >= vfi);
    EXPECT_TRUE(vfi == vfi);
    EXPECT_TRUE(vfi >= vfi);
    EXPECT_TRUE(vfi < vfi8);
    EXPECT_TRUE(vfi <= vfi8);

    // Offset
    EXPECT_EQ(9, vfi[2]);

    // Change the value
    *vfi = 100;
    EXPECT_EQ(100, *vfi);
}

//---------------------------------------------------------------------------//
// const_VF_Iterator TESTS
//---------------------------------------------------------------------------//

TEST_F(View_Field_IteratorTest, const_stride_1_test)
{
    // Default constructor
    cVFI def_vfi;
    EXPECT_EQ(1, def_vfi.stride());
    EXPECT_FALSE(def_vfi.is_valid());

    // Constructor
    cVFI vfi(vi.data(), 1);
    EXPECT_EQ(1, vfi.stride());
    EXPECT_EQ(1, *(vfi.get_pointer()));
    EXPECT_TRUE(vfi.is_valid());

    // Copy constructor
    cVFI vfi2(vfi);
    EXPECT_EQ(1, vfi2.stride());
    EXPECT_EQ(1, *(vfi2.get_pointer()));

    // VF_Iterator -> const_VF_Iterator constructor
    VFI tmp(vi.data(), 1);
    cVFI vfi2_5(tmp);
    EXPECT_EQ(1, *vfi2_5);

    // Copy operator
    cVFI vfi3;
    vfi3 = vfi2;
    EXPECT_EQ(1, vfi3.stride());
    EXPECT_EQ(1, *(vfi3.get_pointer()));

    // Dereference
    EXPECT_EQ(1, *vfi);
    // Prefix increment
    ++vfi;
    EXPECT_EQ(2, *vfi);
    // Prefix decrement
    --vfi;
    EXPECT_EQ(1, *vfi);

    // Postfix increment
    cVFI vfi4 = vfi++;
    EXPECT_EQ(1, *vfi4);
    EXPECT_EQ(2, *vfi);
    // Postfix decrement
    cVFI vfi5 = vfi--;
    EXPECT_EQ(2, *vfi5);
    EXPECT_EQ(1, *vfi);

    // Addition operator
    cVFI vfi6 = vfi+2;
    EXPECT_EQ(3, *vfi6);
    // Subtraction operator
    cVFI vfi7 = vfi6-2;
    EXPECT_EQ(1, *vfi7);
    EXPECT_EQ(2, vfi6 - vfi7);

    // Compound addition
    vfi += 2;
    EXPECT_EQ(3, *vfi);
    // Compound subtraction
    vfi -= 2;
    EXPECT_EQ(1, *vfi);

    // Comparators
    cVFI vfi8 = vfi + 1;
    EXPECT_TRUE(vfi8 > vfi);
    EXPECT_TRUE(vfi8 >= vfi);
    EXPECT_TRUE(vfi == vfi);
    EXPECT_TRUE(vfi >= vfi);
    EXPECT_TRUE(vfi < vfi8);
    EXPECT_TRUE(vfi <= vfi8);

    // Offset
    EXPECT_EQ(3, vfi[2]);
}

//---------------------------------------------------------------------------//

TEST_F(View_Field_IteratorTest, const_stride_4_test)
{
    // Default constructor
    cVFI def_vfi;
    EXPECT_EQ(1, def_vfi.stride());

    // Constructor
    cVFI vfi(vi.data(), 4);
    EXPECT_EQ(4, vfi.stride());
    EXPECT_EQ(1, *(vfi.get_pointer()));

    // Copy constructor
    cVFI vfi2(vfi);
    EXPECT_EQ(4, vfi2.stride());
    EXPECT_EQ(1, *(vfi2.get_pointer()));

    // Copy operator
    cVFI vfi3;
    vfi3 = vfi2;
    EXPECT_EQ(4, vfi3.stride());
    EXPECT_EQ(1, *(vfi3.get_pointer()));

    // Dereference
    EXPECT_EQ(1, *vfi);
    // Prefix increment
    ++vfi;
    EXPECT_EQ(5, *vfi);
    // Prefix decrement
    --vfi;
    EXPECT_EQ(1, *vfi);

    // Postfix increment
    cVFI vfi4 = vfi++;
    EXPECT_EQ(1, *vfi4);
    EXPECT_EQ(5, *vfi);
    // Postfix decrement
    cVFI vfi5 = vfi--;
    EXPECT_EQ(5, *vfi5);
    EXPECT_EQ(1, *vfi);

    // Addition operator
    cVFI vfi6 = vfi+2;
    EXPECT_EQ(9, *vfi6);
    // Subtraction operator
    cVFI vfi7 = vfi6-2;
    EXPECT_EQ(1, *vfi7);
    EXPECT_EQ(2, vfi6 - vfi7);

    // Compound addition
    vfi += 2;
    EXPECT_EQ(9, *vfi);
    // Compound subtraction
    vfi -= 2;
    EXPECT_EQ(1, *vfi);

    // Comparators
    cVFI vfi8 = vfi + 1;
    EXPECT_TRUE(vfi8 > vfi);
    EXPECT_TRUE(vfi8 >= vfi);
    EXPECT_TRUE(vfi == vfi);
    EXPECT_TRUE(vfi >= vfi);
    EXPECT_TRUE(vfi < vfi8);
    EXPECT_TRUE(vfi <= vfi8);

    // Offset
    EXPECT_EQ(9, vfi[2]);
}

//---------------------------------------------------------------------------//
//                    end of tstView_Field_Iterator.cc
//---------------------------------------------------------------------------//
