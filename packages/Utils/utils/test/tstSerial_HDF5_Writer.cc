//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/test/tstSerial_HDF5_Writer.cc
 * \author Thomas M. Evans
 * \date   Thu Nov 07 22:29:22 2013
 * \brief  Unit test of HDF5_IO.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../Serial_HDF5_Writer.hh"

//---------------------------------------------------------------------------//

class HDF5_IO_Test : public testing::Test
{
  protected:
    typedef profugus::Serial_HDF5_Writer IO_t;
    typedef IO_t::Vec_Int                Vec_Int;
    typedef IO_t::Vec_Dbl                Vec_Dbl;
    typedef IO_t::Vec_Char               Vec_Char;
    typedef IO_t::std_string             std_string;

    void SetUp()
    {
        nodes = profugus::nodes();
        node  = profugus::node();
    }

    void add_data()
    {
        double     x = 1.1;
        int        y = 3;
        std_string s = "hello";
        char       c = 'a';

        int    ints[4]  = {1, 3, 5, 7};
        double dbls[3]  = {1.1, 2.4, 10.1};
        char   chars[2] = {'k', 'm'};

        Vec_Int  i(ints, ints + 4);
        Vec_Dbl  d(dbls, dbls + 3);
        Vec_Char cs(chars, chars + 2);

        io.write("a_double", x);
        io.write("an_int", y);
        io.write("a_string", s);
        io.write("a_char", c);

        io.write("a_vec_int", i);
        io.write("a_vec_dbl", d);
        io.write("a_vec_char", cs);
    }

    void open(const std_string &filenm)
    {
        file = H5Fopen(filenm.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    }

    void close()
    {
        H5Fclose(file);
    }

    void read(const std_string &name,
              double            &value)
    {
        H5LTread_dataset_double(file, name.c_str(), &value);
    }

    void read(const std_string &name,
              int               &value)
    {
        H5LTread_dataset_int(file, name.c_str(), &value);
    }

    void read(const std_string &name,
              std_string       &value)
    {
        // hdf5 parameters
        int         rank  = 0;
        size_t      bytes = 0;
        H5T_class_t dt_class;

        H5LTget_dataset_ndims(file, name.c_str(), &rank);
        hsize_t ndims[1] = {0};
        H5LTget_dataset_info(file, name.c_str(), ndims, &dt_class, &bytes);
        std::vector<char> b(bytes);
        H5LTread_dataset_string(file, name.c_str(), &b[0]);
        value = std_string(&b[0]);
    }

    void read(const std_string &name,
              char             &value)
    {
        H5LTread_dataset_char(file, name.c_str(), &value);
    }

    void read(const std_string &name,
              Vec_Dbl           &value)
    {
        // hdf5 parameters
        int         rank  = 0;
        size_t      bytes = 0;
        H5T_class_t dt_class;

        H5LTget_dataset_ndims(file, name.c_str(), &rank);
        hsize_t ndims[1] = {0};
        H5LTget_dataset_info(file, name.c_str(), ndims, &dt_class, &bytes);
        value.resize(ndims[0]);
        H5LTread_dataset_double(file, name.c_str(), &value[0]);
    }

    void read(const std_string &name,
              Vec_Int           &value)
    {
        // hdf5 parameters
        int         rank  = 0;
        size_t      bytes = 0;
        H5T_class_t dt_class;

        H5LTget_dataset_ndims(file, name.c_str(), &rank);
        hsize_t ndims[1] = {0};
        H5LTget_dataset_info(file, name.c_str(), ndims, &dt_class, &bytes);
        value.resize(ndims[0]);
        H5LTread_dataset_int(file, name.c_str(), &value[0]);
    }

    void read(const std_string &name,
              Vec_Char         &value)
    {
        // hdf5 parameters
        int         rank  = 0;
        size_t      bytes = 0;
        H5T_class_t dt_class;

        H5LTget_dataset_ndims(file, name.c_str(), &rank);
        hsize_t ndims[1] = {0};
        H5LTget_dataset_info(file, name.c_str(), ndims, &dt_class, &bytes);
        value.resize(ndims[0]);
        H5LTread_dataset_char(file, name.c_str(), &value[0]);
    }

  protected:
    IO_t io;

    hid_t file;

    int node, nodes;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(HDF5_IO_Test, Open_Close)
{
    io.open("hdf5_oc.h5");
    io.close();
    // double close
    io.close();
}

//---------------------------------------------------------------------------//

TEST_F(HDF5_IO_Test, Write)
{
    // write data
    {
        io.open("hdf5_rw.h5");

        add_data();

        io.close();
    }

    // read data
    if (node == 0)
    {
        open("hdf5_rw.h5");

        int x    = 0;
        double y = 0.0;
        char   c;
        std_string s;

        Vec_Int vi(12);
        Vec_Dbl vd;
        Vec_Char vc;

        read("an_int", x);
        read("a_double", y);
        read("a_char", c);

        read("a_vec_int", vi);
        read("a_vec_dbl", vd);
        read("a_vec_char", vc);

        read("a_string", s);

        EXPECT_EQ(1.1, y);
        EXPECT_EQ(3, x);
        EXPECT_EQ(std_string("hello"), s);
        EXPECT_EQ('a', c);

        Vec_Int ri(4);
        Vec_Dbl rd(3);
        Vec_Char rc(2);
        ri[0] = 1; ri[1] = 3, ri[2] = 5; ri[3] = 7;
        rd[0] = 1.1; rd[1] = 2.4; rd[2] = 10.1;
        rc[0] = 'k'; rc[1] = 'm';

        EXPECT_VEC_EQ(ri, vi);
        EXPECT_VEC_EQ(rd, vd);
        EXPECT_EQ(rc[0], vc[0]);
        EXPECT_EQ(rc[1], vc[1]);

        close();
    }
}

//---------------------------------------------------------------------------//

TEST_F(HDF5_IO_Test, Query)
{
    io.open("hdf5_q.h5");

    add_data();

    if (node == 0)
    {
        EXPECT_TRUE(io.query("an_int"));
        EXPECT_TRUE(io.query("a_double"));
        EXPECT_TRUE(io.query("a_string"));
        EXPECT_TRUE(io.query("a_char"));
        EXPECT_TRUE(io.query("a_vec_int"));
        EXPECT_TRUE(io.query("a_vec_dbl"));
        EXPECT_TRUE(io.query("a_vec_char"));

        EXPECT_FALSE(io.query("a_foo"));
    }
    else
    {
        EXPECT_FALSE(io.query("an_int"));
        EXPECT_FALSE(io.query("a_double"));
        EXPECT_FALSE(io.query("a_string"));
        EXPECT_FALSE(io.query("a_char"));
        EXPECT_FALSE(io.query("a_vec_int"));
        EXPECT_FALSE(io.query("a_vec_dbl"));
        EXPECT_FALSE(io.query("a_vec_char"));
    }

    io.close();
}

//---------------------------------------------------------------------------//

TEST_F(HDF5_IO_Test, Append)
{
    // make basic file
    {
        io.open("hdf5_a.h5");

        add_data();

        io.close();
    }

    // now append
    {
        io.open("hdf5_a.h5", IO_t::APPEND);

        io.write("another_int", 5);

        io.close();
    }

    // read
    if (node == 0)
    {
        open("hdf5_a.h5");

        int x = 0;

        read("another_int", x);

        EXPECT_EQ(5, x);

        close();
    }
}

//---------------------------------------------------------------------------//

TEST_F(HDF5_IO_Test, Non_Master_Write)
{
    if (nodes != 2)
        return;

    io.open("hdf5_node_1.h5", IO_t::CLOBBER, 1);

    int x = 10;
    if (node == 1)
    {
        x = 11;
    }
    io.write("x", x);
    io.close();

    profugus::global_barrier();

    // read
    if (node == 0)
    {
        open("hdf5_node_1.h5");
        int x = 0;
        read("x", x);
        EXPECT_EQ(11, x);
        close();
    }
}

//---------------------------------------------------------------------------//

TEST_F(HDF5_IO_Test, Groups)
{
    // add groups
    io.open("hdf5_groups.h5");
    if (node == 0)
    {
        EXPECT_EQ("/", io.abspath());
    }
    else
    {
        EXPECT_EQ("", io.abspath());
    }

    io.begin_group("dir1");
    if (node == 0)
    {
        EXPECT_EQ("/dir1/", io.abspath());
    }
    else
    {
        EXPECT_EQ("", io.abspath());
    }
    add_data();
    io.end_group();
    if (node == 0)
    {
        EXPECT_EQ("/", io.abspath());
    }
    else
    {
        EXPECT_EQ("", io.abspath());
    }
    io.begin_group("dir2");
    io.begin_group("subdir");
    if (node == 0)
    {
        EXPECT_EQ("/dir2/subdir/", io.abspath());
    }
    else
    {
        EXPECT_EQ("", io.abspath());
    }
    add_data();
    io.end_group();
    io.end_group();

    io.close();
}

//---------------------------------------------------------------------------//
//                 end of tstSerial_HDF5_Writer.cc
//---------------------------------------------------------------------------//
