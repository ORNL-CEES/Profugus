//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/test/tstPacking_Utils.cc
 * \author Thomas M. Evans
 * \date   Fri Apr 25 15:19:23 2014
 * \brief  Packing_Utils test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Packing_Utils.hh"

#include "gtest/utils_gtest.hh"
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <typeinfo>
#include <cstdlib>

using profugus::Packer;
using profugus::Unpacker;
using profugus::pack_data;
using profugus::unpack_data;

using namespace std;

void do_some_packing(Packer &p,
                     const vector<double> &vd,
                     const vector<int> &vi)
{
    for ( int i = 0; i < vd.size(); ++i )
        p << vd[i];

    for ( int i = 0; i < vi.size(); ++i )
        p << vi[i];
}

//---------------------------------------------------------------------------//
// TEST HARNESS
//---------------------------------------------------------------------------//

class PackingUtilsTest : public ::testing::Test
{
  protected:
};

//---------------------------------------------------------------------------//

TEST_F(PackingUtilsTest, compute_buffer_size)
{
    // make data

    int num_vd = 5;
    vector<double> vd(num_vd, 2.3432);
    vd[3] = 22.4;

    int num_vi = 3;
    vector<int> vi(num_vi, 6);
    vi[0] = 7;
    vi[1] = 22;

    unsigned int total_size = num_vi * sizeof(int) + num_vd * sizeof(double);

    Packer p;

    // Compute the required buffer size.

    p.compute_buffer_size_mode();
    do_some_packing(p, vd, vi); // computes the size

    EXPECT_EQ(p.size(), total_size);

    vector<char> buffer(p.size());

    // Pack into buffer.

    p.set_buffer(p.size(), &buffer[0]);
    do_some_packing(p, vd, vi); // actually does the packing

    // Unpack

    Unpacker u;

    u.set_buffer(p.size(), &buffer[0]);

    for ( int i = 0; i < vd.size(); ++i )
    {
        double d;
        u >> d;
        if ( ! soft_equiv(d, vd[i]) ) ADD_FAILURE();
    }

    for ( int i = 0; i < vi.size(); ++i )
    {
        int j;
        u >> j;
        EXPECT_EQ(vi[i], j);
    }
}

//---------------------------------------------------------------------------//

void packing_test()
{
    // make some data
    double x = 102.45;
    double y = 203.89;
    double z = 203.88;

    int ix   = 10;
    int iy   = 11;
    int iz   = 12;

    // make 2 buffers for data
    int s1   = 2 * sizeof(double) + 2 * sizeof(int);
    char *b1 = new char[s1];
    int s2   = sizeof(double) + sizeof(int);
    char *b2 = new char[s2];

    // pack the data
    {
        Packer p;

        p.set_buffer(s1, b1);
        p << x << ix;
        p.pack(y);
        p.pack(iy);

        EXPECT_EQ(b1+s1, p.get_ptr());

        // try catching a failure
#ifdef ENSURE_ON
        EXPECT_THROW(p << iz, profugus::assertion);
#endif
        p.set_buffer(s2, b2);
        p << iz << z;

        EXPECT_EQ(b2+s2, p.get_ptr());
    }

    // unpack the data
    {
        Unpacker u;

        double   d = 0;
        int      i = 0;

        u.set_buffer(s1, b1);
        u >> d >> i;
        EXPECT_EQ(102.45, d);
        EXPECT_EQ(10, i);

        u.unpack(d);
        u.unpack(i);
        EXPECT_EQ(203.89, d);
        EXPECT_EQ(11, i);

        EXPECT_EQ(s1+b1, u.get_ptr());

#ifdef ENSURE_ON
        EXPECT_THROW(u.unpack(i), profugus::assertion);
#endif

        u.set_buffer(s2, b2);
        u >> i >> d;
        EXPECT_EQ(12, i);
        EXPECT_EQ(203.88, d);

        EXPECT_EQ(s2+b2, u.get_ptr());
    }

    delete [] b1;
    delete [] b2;

    // try packing a vector and char array
    double r = 0;
    srand(125);

    vector<double> vx(100, 0.0);
    vector<double> ref(100, 0.0);

    char c[4] = {'c','h','a','r'};

    for (int i = 0; i < vx.size(); i++)
    {
        r      = rand();
        vx[i]  = r / RAND_MAX;
        ref[i] = vx[i];
    }

    int size     = 100 * sizeof(double) + 4;
    char *buffer = new char[size];

    // pack
    {
        Packer p;
        p.set_buffer(size, buffer);

        for (int i = 0; i < vx.size(); i++)
            p << vx[i];

        for (int i = 0; i < 4; i++)
            p << c[i];

        EXPECT_EQ(buffer+size, p.get_ptr());
    }

    // unpack
    {
        char cc[4];
        vector<double> x(100, 0.0);

        Unpacker u;
        u.set_buffer(size, buffer);

        for (int i = 0; i < x.size(); i++)
            u >> x[i];

        for (int i = 0; i < 4; i++)
            u >> cc[i];

        EXPECT_EQ(buffer+size, u.get_ptr());

        for (int i = 0; i < x.size(); i++)
            EXPECT_EQ(ref[i], x[i]);

        EXPECT_EQ('c', c[0]);
        EXPECT_EQ('h', c[1]);
        EXPECT_EQ('a', c[2]);
        EXPECT_EQ('r', c[3]);
    }

    delete [] buffer;
}

//---------------------------------------------------------------------------//
TEST_F(PackingUtilsTest, std_string)
{
    vector<char> pack_string;

    {
        // make a string
        string hw("Hello World");

        // make a packer
        Packer packer;

        // make a char to write the string into
        pack_string.resize(hw.size() + 1 * sizeof(int));

        packer.set_buffer(pack_string.size(), &pack_string[0]);

        packer << static_cast<int>(hw.size());

        // pack it
        for (string::const_iterator it = hw.begin(); it != hw.end(); it++)
            packer << *it;

        EXPECT_EQ(&pack_string[0] + pack_string.size(), packer.get_ptr());
        EXPECT_EQ(packer.begin()  + pack_string.size(), packer.get_ptr());
    }

    // now unpack it
    Unpacker unpacker;
    unpacker.set_buffer(pack_string.size(), &pack_string[0]);

    // unpack the size of the string
    int size;
    unpacker >> size;

    string nhw;
    nhw.resize(size);

    // unpack the string
    for (string::iterator it = nhw.begin(); it != nhw.end(); it++)
        unpacker >> *it;

    EXPECT_EQ(&pack_string[0] + pack_string.size(), unpacker.get_ptr());

    // test the unpacked string
    // make a string
    string hw("Hello World");

    if (hw == nhw)
    {
        ostringstream message;
        message << "Unpacked string " << nhw << " that matches original "
                << "string " << hw;
        SUCCEED() << message.str();
    }
    else
    {
        ostringstream message;
        message << "Failed to unpack string " << hw << " correctly. Instead "
                << "unpacked " << nhw;
        FAIL() << message.str();
    }
}

//---------------------------------------------------------------------------//
TEST_F(PackingUtilsTest, packing_functions)
{
    // make a packing container
    vector<char> total_packed;

    vector<double> x_ref;
    string         y_ref("Oak Ridge is a fun place!");

    vector<double> x_new;
    string         y_new;

    // pack some data
    {
        vector<char> packed_int(sizeof(int));
        vector<char> packed_vector;
        vector<char> packed_string;

        vector<double> x(5);
        string         y("Oak Ridge is a fun place!");

        x_ref.resize(5);

        for (int i = 0; i < 5; i++)
        {
            x[i]     = 100.0 * (static_cast<double>(i) + x[i]) + 2.5;
            x_ref[i] = x[i];
        }

        pack_data(x, packed_vector);
        pack_data(y, packed_string);

        EXPECT_EQ(5 * sizeof(double) + sizeof(int), packed_vector.size());
        EXPECT_EQ(y_ref.size() + sizeof(int), packed_string.size());

        Packer p;
        p.set_buffer(sizeof(int), &packed_int[0]);

        // pack the size of the vector
        p << static_cast<int>(packed_vector.size());

        // push the vector onto the total packed
        total_packed.insert(total_packed.end(),
                            packed_int.begin(), packed_int.end());
        total_packed.insert(total_packed.end(),
                            packed_vector.begin(), packed_vector.end());

        // reset the packer
        p.set_buffer(sizeof(int), &packed_int[0]);

        // pack the size of the string
        p << static_cast<int>(packed_string.size());

        // push the string onto the total packed
        total_packed.insert(total_packed.end(),
                            packed_int.begin(), packed_int.end());
        total_packed.insert(total_packed.end(),
                            packed_string.begin(), packed_string.end());

        if (total_packed.size() != 2*packed_int.size() + packed_vector.size()
            + packed_string.size()) ADD_FAILURE();
    }

    // unpack the data
    {
        int size;
        Unpacker u;
        u.set_buffer(total_packed.size(), &total_packed[0]);

        // unpack the packed vector
        u >> size;
        vector<char> packed_vector(size);

        for (int i = 0; i < size; i++)
            u >> packed_vector[i];

        // unpack the packed string
        u >> size;
        vector<char> packed_string(size);

        for (int i = 0; i < size; i++)
            u >> packed_string[i];

        EXPECT_EQ(&total_packed[0] + total_packed.size(), u.get_ptr());

        unpack_data(x_new, packed_vector);
        unpack_data(y_new, packed_string);
    }

    if (!soft_equiv(x_new.begin(), x_new.end(), x_ref.begin(), x_ref.end()))
        ADD_FAILURE();

    if (y_new == y_ref)
    {
        ostringstream message;
        message << "Correctly unpacked string " << y_new << " from "
                << "packed " << y_ref;
        SUCCEED() << message.str();
    }
    else
    {
        ostringstream message;
        message << "Failed to unpack string " << y_ref << " correctly. Instead "
                << "unpacked " << y_new;
        FAIL() << message.str();
    }
}

//---------------------------------------------------------------------------//
//                 end of tstPacking_Utils.cc
//---------------------------------------------------------------------------//
