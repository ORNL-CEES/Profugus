//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/test/tstParallel_HDF5_Writer.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 24 09:51:23 2014
 * \brief  Parallel_HDF5_Writer unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <algorithm>
#include "../Parallel_HDF5_Writer.hh"

#ifdef H5_HAVE_PARALLEL

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class HDF_IO_Test : public testing::Test
{
  protected:
    typedef profugus::Parallel_HDF5_Writer IO_t;
    typedef IO_t::Vec_Int                  Vec_Int;
    typedef IO_t::Vec_Dbl                  Vec_Dbl;
    typedef IO_t::std_string               std_string;
    typedef IO_t::Decomp                   Decomp;
    typedef IO_t::Vec_Hsize                Vec_Hsize;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        nodes = profugus::nodes();
        node  = profugus::node();
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
              double           &value)
    {
        H5LTread_dataset_double(file, name.c_str(), &value);
    }

    void read(const std_string &name,
              int              &value)
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
              Vec_Int          &value)
    {
        // hdf5 parameters
        int         rank  = 0;
        size_t      bytes = 0;
        H5T_class_t dt_class;

        H5LTget_dataset_ndims(file, name.c_str(), &rank);
        Vec_Hsize ndims(rank);
        H5LTget_dataset_info(file, name.c_str(), &ndims[0], &dt_class, &bytes);
        int size = 1;
        for (int n = 0; n < rank; ++n)
        {
            size *= ndims[n];
        }
        value.resize(size);
        H5LTread_dataset_int(file, name.c_str(), &value[0]);
    }

    void read(const std_string &name,
              Vec_Dbl          &value)
    {
        // hdf5 parameters
        int         rank  = 0;
        size_t      bytes = 0;
        H5T_class_t dt_class;

        H5LTget_dataset_ndims(file, name.c_str(), &rank);
        Vec_Hsize ndims(rank);
        H5LTget_dataset_info(file, name.c_str(), &ndims[0], &dt_class, &bytes);
        int size = 1;
        for (int n = 0; n < rank; ++n)
        {
            size *= ndims[n];
        }
        value.resize(size);
        H5LTread_dataset_double(file, name.c_str(), &value[0]);
    }

  protected:
    // >>> Data that get re-initialized between tests

    int node, nodes;

    hid_t file;
};

//---------------------------------------------------------------------------//

class Parallel_Field_IO : public HDF_IO_Test
{
  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        nodes = profugus::nodes();
        node  = profugus::node();
    }

    // Build the spatial decomposition.
    void decompose(int Nx,
                   int Ny,
                   int Nz,
                   int Px,
                   int Py,
                   int Pz)
    {
        global[0]  = Nx;
        global[1]  = Ny;
        global[2]  = Nz;
        num_global = Nx * Ny * Nz;

        EXPECT_EQ(0, Nx % Px);
        EXPECT_EQ(0, Ny % Py);
        EXPECT_EQ(0, Nz % Pz);

        local[0]  = Nx / Px;
        local[1]  = Ny / Py;
        local[2]  = Nz / Pz;
        num_local = local[0] * local[1] * local[2];

        // Spatial decomposition index
        //   idx = i + Nx * j + Nx * Ny * k

        // block id
        int bidx[] = {0, 0, 0};

        bidx[2] = node / (Py * Px);
        bidx[1] = (node - Py*Px*bidx[2]) / Px;
        bidx[0] = node - Px*bidx[1] - Py*Px*bidx[2];

        offset[0] = local[0] * bidx[0];
        offset[1] = local[1] * bidx[1];
        offset[2] = local[2] * bidx[2];

        EXPECT_EQ(node, bidx[0] + Px * (bidx[1] + Py * (bidx[2])));
    }

    // Local-to-global mapping
    int l2l(int i, int j, int k)
    {
        return i + local[0] * (j + local[1] * (k));
    }

    // Local-to-global mapping
    int l2g(int i, int j, int k)
    {
        return (i + offset[0]) + global[0] * (
            (j + offset[1]) + global[1] * (k + offset[2]));
    }

  protected:

    int local[3];
    int global[3];
    int offset[3];

    int num_global, num_local;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(HDF_IO_Test, No_Write)
{
    IO_t writer;
    writer.open("no_write.h5");

    writer.close();
}

//---------------------------------------------------------------------------//

TEST_F(HDF_IO_Test, Decomposition)
{
    Decomp d(5);

    EXPECT_EQ(5, d.ndims);
    EXPECT_EQ(5, d.global.size());
    EXPECT_EQ(5, d.local.size());
    EXPECT_EQ(5, d.offset.size());
    EXPECT_EQ(IO_t::COLUMN_MAJOR, d.order);

    for (int i = 0; i < 5; ++i)
    {
        EXPECT_EQ(0, d.global[i]);
        EXPECT_EQ(0, d.local[i]);
        EXPECT_EQ(0, d.offset[i]);
    }

    // flip data
    d.global[0] = 1;
    d.global[1] = 6;
    d.global[2] = 4;
    d.global[3] = 2;
    d.global[4] = 11;

    IO_t::Vec_Hsize g(5);
    std::copy(d.global.rbegin(), d.global.rend(), g.begin());
    EXPECT_EQ(11, g[0]);
    EXPECT_EQ(2,  g[1]);
    EXPECT_EQ(4,  g[2]);
    EXPECT_EQ(6,  g[3]);
    EXPECT_EQ(1,  g[4]);
}

//---------------------------------------------------------------------------//
/*
 * This test assumes one unknown per mesh cell.  For example, integrated
 * scalar flux.  The field is ordered COLUMN-MAJOR, ie.
 *    n = i + Nx * (j + Ny * (k))
 */
TEST_F(Parallel_Field_IO, OneUnknown_4Blocks)
{
    if (nodes != 4)
        return;

    // decompose a 6x4x2 mesh into 4 blocks
    decompose(6, 4, 2, 2, 2, 1);

    IO_t writer;
    writer.open("one_48B.h5");

    // make a global reference field
    Vec_Dbl ref(48, 0.0);
    for (int n = 0; n < 48; ++n)
        ref[n] = n + 1;

    // make the decomposition (the decomposition is over (i,j,k))
    Decomp d(3);
    std::copy(local,  local + 3,  d.local.begin());
    std::copy(global, global + 3, d.global.begin());
    std::copy(offset, offset + 3, d.offset.begin());

    // write using default COLUMN-MAJOR ORDER
    {
        // write local data
        Vec_Dbl data(num_local, 0);
        for (int k = 0; k < local[2]; ++k)
        {
            for (int j = 0; j < local[1]; ++j)
            {
                for (int i = 0; i < local[0]; ++i)
                {
                    data[l2l(i,j,k)] = ref[l2g(i,j,k)];
                }
            }
        }

        writer.write("mesh_data_CM", d, &data[0]);
    }

    // write using ROW-MAJOR ORDER (we reorder the data to be ROW-MAJOR)
    {
        d.order = IO_t::ROW_MAJOR;

        // write local data
        Vec_Dbl data(num_local, 0);
        for (int k = 0; k < local[2]; ++k)
        {
            for (int j = 0; j < local[1]; ++j)
            {
                for (int i = 0; i < local[0]; ++i)
                {
                    int l = k + local[2] * (j + local[1] * (i));
                    int g = (k+offset[2]) + global[2] * (
                        (j+offset[1]) + global[1] * (i+offset[0]));
                    data[l] = ref[g] + 100.0;
                }
            }
        }

        writer.write("mesh_data_RM", d, &data[0]);
    }

    writer.close();

    profugus::global_barrier();

    // read data
    if (node == 0)
    {
        open("one_48B.h5");

        Vec_Dbl d1, d2;
        read("mesh_data_RM", d2);
        read("mesh_data_CM", d1);

        EXPECT_EQ(48, d1.size());
        EXPECT_EQ(48, d2.size());

        EXPECT_VEC_EQ(ref, d1);

        for (int n = 0; n < 48; ++n)
        {
            EXPECT_EQ(ref[n] + 100.0, d2[n]);
        }

        close();
    }
}

//---------------------------------------------------------------------------//
/*
 * This test assumes 4 unknowns per mesh cell and 3 groups.  The field is
 * ordered COLUMN-MAJOR, ie.
 *   n = u + Nu * (i + Nx * (j + Ny * (k + Nz * (g))))
 *
 * However, the decomposition is only over space.
 */
TEST_F(Parallel_Field_IO, FourUnknowns_4Blocks)
{
    if (nodes != 4)
        return;

    // decompose a 6x4x2 mesh into 4 blocks
    decompose(6, 4, 2, 2, 2, 1);

    IO_t writer;
    writer.open("four_48B.h5");

    // add a group to the file
    writer.begin_group("data");

    // make a global reference field
    Vec_Dbl ref(576, 0.0);
    for (int n = 0; n < 576; ++n)
        ref[n] = n + 1.0;

    // make the decomposition (the decomposition is over (i,j,k)), but the
    // data is tiled over 5 variables
    Decomp d(5);
    d.local[0]  = 4; // 4 unknowns per cell
    d.local[4]  = 3; // 3 groups
    d.global[0] = 4; // 4 unknowns per cell
    d.global[4] = 3; // 3 groups
    d.offset[0] = 0; // unknowns are not decomposed
    d.offset[4] = 0; // groups are not decomposed
    std::copy(local,  local + 3,  d.local.begin() + 1);
    std::copy(global, global + 3, d.global.begin() + 1);
    std::copy(offset, offset + 3, d.offset.begin() + 1);

    // write using default COLUMN-MAJOR ORDER
    {
        // write local data
        Vec_Dbl data(num_local * 12, 0);
        for (int g = 0; g < 3; ++g)
        {
            for (int k = 0; k < local[2]; ++k)
            {
                for (int j = 0; j < local[1]; ++j)
                {
                    for (int i = 0; i < local[0]; ++i)
                    {
                        for (int n = 0; n < 4; ++n)
                        {
                            int ln = n + 4 * l2l(i, j, k) + 4 * num_local * g;
                            int gn = n + 4 * l2g(i, j, k) + 4 * num_global * g;
                            data[ln] = ref[gn];
                        }
                    }
                }
            }
        }

        writer.write("mesh_data", d, &data[0]);
    }

    writer.end_group();

    writer.close();

    profugus::global_barrier();

    // read data
    if (node == 0)
    {
        open("four_48B.h5");

        Vec_Dbl d;
        read("/data/mesh_data", d);

        EXPECT_EQ(576, d.size());

        EXPECT_VEC_EQ(ref, d);

        close();
    }
}

//---------------------------------------------------------------------------//

#else

TEST(HDF5_IO, Parallel_HDF5_Unavailable)
{
}

#endif // H5_HAVE_PARALLEL

//---------------------------------------------------------------------------//
//                 end of tstParallel_HDF5_Writer.cc
//---------------------------------------------------------------------------//
