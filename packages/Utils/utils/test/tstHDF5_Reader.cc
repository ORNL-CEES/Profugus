//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/test/tstHDF5_Reader.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 24 09:51:32 2014
 * \brief  HDF5_Reader unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../HDF5_Reader.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class HDF5_IO_Test : public testing::Test
{
  protected:
    typedef denovo::HDF5_Reader IO_t;
    typedef IO_t::Vec_Int       Vec_Int;
    typedef IO_t::Vec_Dbl       Vec_Dbl;
    typedef IO_t::Vec_Char      Vec_Char;
    typedef IO_t::std_string    std_string;
    typedef IO_t::Decomp        Decomp;
    typedef IO_t::Vec_Hsize     Vec_Hsize;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        nodes = profugus::nodes();
        node  = profugus::node();

        // write some data fields
        if (node == 0)
        {
            hid_t file = H5Fcreate("hdf5_reader.h5", H5F_ACC_TRUNC,
                                   H5P_DEFAULT, H5P_DEFAULT);
            hid_t group = H5Gcreate(file, "scalar", H5P_DEFAULT, H5P_DEFAULT,
                                    H5P_DEFAULT);

            // write some ints
            int x = 10, y = 11, z = 12;
            hsize_t dims[] = {1};

            H5LTmake_dataset_int(group, "x", 1, dims, &x);
            H5LTmake_dataset_int(group, "y", 1, dims, &y);
            H5LTmake_dataset_int(group, "z", 1, dims, &z);

            // write some chars
            char c1 = 'a', c2 = 'b', c3 = 'z';

            H5LTmake_dataset_char(group, "c1", 1, dims, &c1);
            H5LTmake_dataset_char(group, "c2", 1, dims, &c2);
            H5LTmake_dataset_char(group, "c3", 1, dims, &c3);

            // write a field
            hsize_t dims3[] = {3, 4, 6};
            hsize_t dims5[] = {5, 3, 4, 6, 4};
            Vec_Dbl data3(3*4*6, 0.0);
            Vec_Dbl data5(5*3*4*6*4, 0.0);
            for (int g = 0; g < 5; ++g)
            {
                for (int k = 0; k < 3; ++k)
                {
                    for (int j = 0; j < 4; ++j)
                    {
                        for (int i = 0; i < 6; ++i)
                        {
                            int c    = i + 6 * (j + 4 * (k));
                            data3[c] = c + 10.0;

                            for (int m = 0; m < 4; ++m)
                            {
                                int n    = m + 4 * (c + 72 * (g));
                                data5[n] = n + 100.0;
                            }
                        }
                    }
                }
            }

            H5LTmake_dataset_double(file, "data3", 3, &dims3[0], &data3[0]);
            H5LTmake_dataset_double(file, "data5", 5, &dims5[0], &data5[0]);

            // write a char field
            hsize_t dims1[] = {5};
            hsize_t dims2[] = {2, 3};
            Vec_Char char1(5);
            Vec_Char char2(2*3);
            char1[0] = 'a'; char1[1] = 'c'; char1[2] = 'e'; char1[3] = 'g';
            char1[4] = 'i';
            char2[0] = 'k'; char2[1] = 'm'; char2[2] = 'o';
            char2[3] = 'q'; char2[4] = 's'; char2[5] = 'u';

            H5LTmake_dataset_char(file, "char1", 1, &dims1[0], &char1[0]);
            H5LTmake_dataset_char(file, "char2", 2, &dims2[0], &char2[0]);

            H5Gclose(group);
            H5Fclose(file);
        }

        profugus::global_barrier();
        reader.open("hdf5_reader.h5", node);
    }

    void TearDown()
    {
        reader.close();
        profugus::global_barrier();
    }

  protected:
    // >>> Data that get re-initialized between tests

    int node, nodes;

    IO_t reader;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(HDF5_IO_Test, read_simple)
{
    EXPECT_FALSE(reader.query("x"));
    EXPECT_FALSE(reader.query("y"));
    EXPECT_FALSE(reader.query("z"));
    EXPECT_FALSE(reader.query("c1"));
    EXPECT_FALSE(reader.query("c2"));
    EXPECT_FALSE(reader.query("c3"));

    reader.begin_group("scalar");

    EXPECT_TRUE(reader.query("x"));
    EXPECT_TRUE(reader.query("y"));
    EXPECT_TRUE(reader.query("z"));

    int x = 0, y = 0, z = 0;
    char c1, c2, c3;

    reader.read("x", x);
    reader.read("y", y);
    reader.read("z", z);
    reader.read("c1", c1);
    reader.read("c2", c2);
    reader.read("c3", c3);

    EXPECT_EQ(10, x);
    EXPECT_EQ(11, y);
    EXPECT_EQ(12, z);
    EXPECT_EQ('a', c1);
    EXPECT_EQ('b', c2);
    EXPECT_EQ('z', c3);

    reader.end_group();
}

//---------------------------------------------------------------------------//

TEST_F(HDF5_IO_Test, fields)
{
    Vec_Dbl d3, d5;
    reader.read("data3", d3);
    reader.read("data5", d5);

    EXPECT_EQ(72, d3.size());
    EXPECT_EQ(1440, d5.size());

    for (int n = 0; n < 72; ++n)
    {
        EXPECT_EQ(n + 10.0, d3[n]);
    }

    for (int n = 0; n < 1440; ++n)
    {
        EXPECT_EQ(n + 100.0, d5[n]);
    }

    Vec_Char vc1, vc2;
    reader.read("char1", vc1);
    reader.read("char2", vc2);

    Vec_Char char1(5);
    Vec_Char char2(2*3);
    char1[0] = 'a'; char1[1] = 'c'; char1[2] = 'e'; char1[3] = 'g';
    char1[4] = 'i';
    char2[0]= 'k'; char2[1] = 'm'; char2[2] = 'o';
    char2[3] = 'q'; char2[4] = 's'; char2[5] = 'u';

    EXPECT_EQ(5, vc1.size());
    EXPECT_EQ(6, vc2.size());

    for (int n = 0; n < 5; ++n)
    {
        EXPECT_EQ(char1[n], vc1[n]);
    }
    for (int n = 0; n < 6; ++n)
    {
        EXPECT_EQ(char2[n], vc2[n]);
    }
}

//---------------------------------------------------------------------------//

TEST_F(HDF5_IO_Test, decomposition)
{
    Decomp d3, d5;
    reader.get_decomposition("data3", d3);
    reader.get_decomposition("data5", d5);

    EXPECT_EQ(3, d3.ndims);
    EXPECT_EQ(6, d3.global[0]);
    EXPECT_EQ(4, d3.global[1]);
    EXPECT_EQ(3, d3.global[2]);

    EXPECT_EQ(5, d5.ndims);
    EXPECT_EQ(4, d5.global[0]);
    EXPECT_EQ(6, d5.global[1]);
    EXPECT_EQ(4, d5.global[2]);
    EXPECT_EQ(3, d5.global[3]);
    EXPECT_EQ(5, d5.global[4]);

    EXPECT_EQ(IO_t::COLUMN_MAJOR, d3.order);
    EXPECT_EQ(IO_t::COLUMN_MAJOR, d5.order);

    for (int n = 0; n < 3; ++n)
    {
        EXPECT_EQ(0, d3.local[n]);
        EXPECT_EQ(0, d3.offset[n]);
    }

    for (int n = 0; n < 5; ++n)
    {
        EXPECT_EQ(0, d5.local[n]);
        EXPECT_EQ(0, d5.offset[n]);
    }
}

//---------------------------------------------------------------------------//

TEST_F(HDF5_IO_Test, chunks)
{
    Decomp d3, d5;
    reader.get_decomposition("data3", d3);
    reader.get_decomposition("data5", d5);

    // get chunks of both fields

    // 3-dim field
    d3.local[0]  = 2;
    d3.local[1]  = 4;
    d3.local[2]  = 3;
    d3.offset[0] = 1;
    d3.offset[1] = 0;
    d3.offset[2] = 0;
    Vec_Dbl data3(24, 0.0);
    reader.read("data3", d3, &data3[0]);

    int ctr = 0;
    for (int k = 0; k < 3; ++k)
    {
        for (int j = 0; j < 4; ++j)
        {
            for (int i = 0; i < 2; ++i)
            {
                int n = (i+1) + 6 * (j + 4 * (k));
                EXPECT_EQ(n + 10.0, data3[ctr]);
                ++ctr;
            }
        }
    }
    EXPECT_EQ(24, ctr);

    // 5-dim field
    d5.local[0]  = 4;
    d5.local[1]  = 2;
    d5.local[2]  = 1;
    d5.local[3]  = 1;
    d5.local[4]  = 5;
    d5.offset[1] = 0;
    d5.offset[2] = 3;
    d5.offset[3] = 1;
    Vec_Dbl data5(40, 0.0);
    reader.read("data5", d5, &data5[0]);

    ctr = 0;
    for (int g = 0; g < 5; ++g)
    {
        for (int k = 0; k < 1; ++k)
        {
            for (int j = 0; j < 1; ++j)
            {
                for (int i = 0; i < 2; ++i)
                {
                    for (int m = 0; m < 4; ++m)
                    {
                        int c = i + 6 * ((j+3) + 4 * (k+1));
                        int n = m + 4 * (c + 72 * (g));
                        EXPECT_EQ(n + 100.0, data5[ctr]);
                        ++ctr;
                    }
                }
            }
        }
    }
    EXPECT_EQ(40, ctr);
}

//---------------------------------------------------------------------------//
//                 end of tstHDF5_Reader.cc
//---------------------------------------------------------------------------//
