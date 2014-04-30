//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/test/tstHDF5.cc
 * \author Thomas M. Evans
 * \date   Thu Dec 12 16:37:32 2013
 * \brief  HDF5 Functionality test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "harness/DBC.hh"
#include "comm/global.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class Block_Data_Simple : public testing::Test
{
  protected:
    typedef std::vector<int> Field;
    typedef std::string      string_t;

  protected:
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();
    }

    // Build the data
    void build(int Nx,
               int Ny,
               int Nz,
               int Px,
               int Py,
               int Pz)
    {
        // dimensions of mesh data -> our 1D data is aligned in COLUMN-MAJOR
        // ordering (FORTRAN), ie. index = i + Nx * (j + Ny * k)
        dims[0] = Nz;
        dims[1] = Ny;
        dims[2] = Nx;
        gcells  = Nx * Ny * Nz;

        Check (Nx % Px == 0);
        Check (Ny % Py == 0);
        Check (Nz % Pz == 0);

        local[2] = Nx / Px;
        local[1] = Ny / Py;
        local[0] = Nz / Pz;
        lcells   = local[0] * local[1] * local[2];

        // num blocks by dim
        int nb[]       = {Pz, Py, Px};
        int num_blocks = Px * Py * Pz;

        // block id
        int block  = node;
        int bidx[] = {0, 0, 0};

        // idx = i + Nbi * j + Nbi * Nbj * k

        bidx[0] = block / (nb[1] * nb[2]);
        bidx[1] = (block - nb[1]*nb[2]*bidx[0]) / nb[2];
        bidx[2] = block - nb[2]*bidx[1] - nb[1]*nb[2]*bidx[0];

        offset[0] = local[0] * bidx[0];
        offset[1] = local[1] * bidx[1];
        offset[2] = local[2] * bidx[2];

        domain.resize(lcells, node);
        data.resize(lcells, 0);

        lids.resize(lcells, 0);
        gids.resize(lcells, 0);
        for (int k = 0; k < local[0]; ++k)
        {
            for (int j = 0; j < local[1]; ++j)
            {
                for (int i = 0; i < local[2]; ++i)
                {
                    int l   = i + local[2] * (j + local[1] * k);
                    lids[l] = l;
                    gids[l] = (i+offset[2]) + Nx * (
                        (j+offset[1]) + Ny * (k+offset[0]));

                    data[l] = 100 + gids[l];
                }
            }
        }
    }

  protected:
    int node, nodes;

    hsize_t dims[3];
    int     gcells;
    hsize_t local[3];
    int     lcells;
    hsize_t offset[3];

    Field lids, gids;
    Field domain;
    Field data;

    hid_t  file_id;
    herr_t status;

    string_t filename;
};

//---------------------------------------------------------------------------//
// Test writing complex multi-D data using simple 1-D approach
// Assume data ordering:
/*
    ind = m + Nm * (i + Nx * (j + Ny * (k + Nz * (g))))
 */
// And decomposition is in (x,y,z)

class Block_Data : public testing::Test
{
  protected:
    typedef std::vector<int> Field;
    typedef std::string      string_t;

  protected:
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();
    }

    // Build the data
    void build(int Nx,
               int Ny,
               int Nz,
               int groups,
               int moments,
               int Px,
               int Py,
               int Pz)
    {
        Nc[0]   = Nx;
        Nc[1]   = Ny;
        Nc[2]   = Nz;
        gcells  = Nx * Ny * Nz;

        Check (Nx % Px == 0);
        Check (Ny % Py == 0);
        Check (Nz % Pz == 0);

        Nl[0]  = Nx / Px;
        Nl[1]  = Ny / Py;
        Nl[2]  = Nz / Pz;
        lcells = Nl[0] * Nl[1] * Nl[2];

        Nm = moments;
        Ng = groups;

        // num blocks by dim
        int nb[]       = {Px, Py, Pz};
        int num_blocks = Px * Py * Pz;

        // block id
        int block  = node;
        int bidx[] = {0, 0, 0};

        // idx = i + Nbi * j + Nbi * Nbj * k

        bidx[2] = block / (nb[1] * nb[0]);
        bidx[1] = (block - nb[1]*nb[0]*bidx[2]) / nb[0];
        bidx[0] = block - nb[0]*bidx[1] - nb[1]*nb[0]*bidx[2];

        coff[0] = Nl[0] * bidx[0];
        coff[1] = Nl[1] * bidx[1];
        coff[2] = Nl[2] * bidx[2];

        // data size
        lsize = lcells * Nm * Ng;
        gsize = gcells * Nm * Ng;

        domain.resize(lsize, node);
        data.resize(lsize, 0);
        lids.resize(lsize, 0);
        gids.resize(lsize, 0);

        for (int g = 0; g < Ng; ++g)
        {
            for (int k = 0; k < Nl[2]; ++k)
            {
                for (int j = 0; j < Nl[1]; ++j)
                {
                    for (int i = 0; i < Nl[0]; ++i)
                    {
                        for (int m = 0; m < Nm; ++m)
                        {
                            int l   = l2l(i, j, k, m, g);
                            lids[l] = l;
                            gids[l] = l2g(i, j, k, m, g);
                            data[l] = 100 + gids[l];
                        }
                    }
                }
            }
        }

        dims[4] = Nm;
        dims[3] = Nc[0];
        dims[2] = Nc[1];
        dims[1] = Nc[2];
        dims[0] = Ng;

        local[4] = Nm;
        local[3] = Nl[0];
        local[2] = Nl[1];
        local[1] = Nl[2];
        local[0] = Ng;

        offset[4] = 0;
        offset[3] = coff[0];
        offset[2] = coff[1];
        offset[1] = coff[2];
        offset[0] = 0;
    }

    int l2l(int i,
            int j,
            int k,
            int m,
            int g)
    {
        return m + Nm * (i + Nl[0] * (j + Nl[1] * (k + Nl[2] * (g))));
    }

    int l2g(int i,
            int j,
            int k,
            int m,
            int g)
    {
        return m + Nm * ((i+coff[0]) + Nc[0] * (
                             (j+coff[1]) + Nc[1] * (
                                 (k+coff[2]) + Nc[2] * (g))));
    }

  protected:

    int Nc[3], Nl[3], Nm, Ng;
    int gcells, lcells;
    int coff[3];
    int lsize, gsize;

    int node, nodes;

    hsize_t dims[5];
    hsize_t offset[5];
    hsize_t local[5];

    Field lids, gids;
    Field domain;
    Field data;

    hid_t  file_id;
    herr_t status;

    string_t filename;
};

//---------------------------------------------------------------------------//

#ifdef H5_HAVE_PARALLEL

class Parallel_HDF5_Simple : public Block_Data_Simple
{
  protected:

    void open(const string_t &basename)
    {
        // make filename
        std::ostringstream f;
        f << basename << "_" << nodes << ".h5";
        filename = f.str();

        // set property list with parallel I/O access
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, profugus::communicator, MPI_INFO_NULL);

        // create a new file collectively and release the file identifier
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            plist_id);
        H5Pclose(plist_id);
    }

    void close()
    {
        H5Fclose(file_id);
    }

    void write(const string_t &dataname,
               const int      *data)
    {
        // create the dataspace for the data
        hid_t filespace = H5Screate_simple(3, dims, NULL);
        hid_t memspace  = H5Screate_simple(3, local, NULL);

        // create the datasets
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 3, local);
        hid_t dset_id = H5Dcreate(file_id, dataname.c_str(), H5T_NATIVE_INT,
                                  filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        H5Pclose(plist_id);
        H5Sclose(filespace);

        // hyperslab parameters
        hsize_t count[]  = {1,1,1}; // 1-block in each dimension
        hsize_t stride[] = {1,1,1}; // stride 1 element at a time

        // select the hyperslace
        filespace = H5Dget_space(dset_id);
        status    = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                        offset, stride, count, local);

        // create property list for collective I/O
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // write the data
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                          plist_id, data);

        // close
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Pclose(plist_id);
    }
};

//---------------------------------------------------------------------------//

class Parallel_HDF5 : public Block_Data
{
  protected:

    void open(const string_t &basename)
    {
        // make filename
        std::ostringstream f;
        f << basename << "_" << nodes << ".h5";
        filename = f.str();

        // set property list with parallel I/O access
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, profugus::communicator, MPI_INFO_NULL);

        // create a new file collectively and release the file identifier
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            plist_id);
        H5Pclose(plist_id);
    }

    void close()
    {
        H5Fclose(file_id);
    }

    void write(const string_t &dataname,
               const int      *data)
    {
        // create the dataspace for the data
        hid_t filespace = H5Screate_simple(5, dims, NULL);
        hid_t memspace  = H5Screate_simple(5, local, NULL);

        // create the datasets
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 5, local);
        hid_t dset_id = H5Dcreate(file_id, dataname.c_str(), H5T_NATIVE_INT,
                                  filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        H5Pclose(plist_id);
        H5Sclose(filespace);

        // hyperslab parameters
        hsize_t count[]  = {1,1,1,1,1}; // 1-block in each dimension
        hsize_t stride[] = {1,1,1,1,1}; // stride 1 element at a time

        // select the hyperslace
        filespace = H5Dget_space(dset_id);
        status    = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                        offset, stride, count, local);

        // create property list for collective I/O
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // write the data
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                          plist_id, data);

        // close
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Pclose(plist_id);
    }

    void simple_write(const string_t &dataname,
                      int             proc,
                      int             size,
                      const int      *data)
    {
        hsize_t dims[]  = {static_cast<hsize_t>(size)};
        hid_t filespace = H5Screate_simple(1, dims, NULL);

        // create the datasets
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, dims);
        hid_t dset_id = H5Dcreate(file_id, dataname.c_str(), H5T_NATIVE_INT,
                                  filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        H5Pclose(plist_id);
        H5Sclose(filespace);

        if (node == proc)
        {
            // hyperslab parameters
            hsize_t count[]  = {1};
            hsize_t stride[] = {1};

            // select the hyperslace
            filespace = H5Dget_space(dset_id);
            status    = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                            offset, stride, count, dims);

            // create property list for collective I/O
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

            // write the data
            status = H5Dwrite(dset_id, H5T_NATIVE_INT,filespace, filespace,
                              plist_id, data);

            H5Sclose(filespace);
            H5Pclose(plist_id);
        }

        // close
        H5Dclose(dset_id);
    }
};

//---------------------------------------------------------------------------//

TEST_F(Parallel_HDF5_Simple, write)
{
    if (nodes != 4)
        return;

    build(6, 4, 2, 2, 2, 1);
    open("write_642");

    write("GIDS", &gids[0]);
    write("LIDS", &lids[0]);
    write("DOMAIN", &domain[0]);
    write("FIELD", &data[0]);

    close();

    // test reading the file
    if (node == 0)
    {
        hid_t rfile = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

        std::vector<int> d(48, 0);

        hsize_t ndims[3]  = {0};
        size_t      bytes = 0;
        H5T_class_t dt_class;
        H5LTget_dataset_info(rfile, "GIDS", ndims, &dt_class, &bytes);

        EXPECT_EQ(2, ndims[0]);
        EXPECT_EQ(4, ndims[1]);
        EXPECT_EQ(6, ndims[2]);

        H5LTread_dataset_int(rfile, "GIDS", &d[0]);

        int c = 0;
        for (int k = 0; k < 2; ++k)
        {
            for (int j = 0; j < 4; ++j)
            {
                for (int i = 0; i < 6; ++i)
                {
                    int gid = i + 6 * (j + 4 * k);
                    EXPECT_EQ(gid, d[c]);
                    ++c;
                }
            }
        }
        EXPECT_EQ(48, c);

        H5LTread_dataset_int(rfile, "LIDS", &d[0]);
        int lids[] = {0,  1,  2, 0,  1,  2,
                      3,  4,  5, 3,  4,  5,
                      0,  1,  2, 0,  1,  2,
                      3,  4,  5, 3,  4,  5,
                      6,  7,  8, 6,  7,  8,
                      9, 10, 11, 9, 10, 11,
                      6,  7,  8, 6,  7,  8,
                      9, 10, 11, 9, 10, 11};

        for (int n = 0; n < 48; ++n)
        {
            EXPECT_EQ(lids[n], d[n]);
        }

        H5Fclose(rfile);
    }
}

//---------------------------------------------------------------------------//

TEST_F(Parallel_HDF5, write)
{
    if (nodes != 4)
        return;

    build(6, 4, 2, 3, 1, 2, 2, 1);
    open("write_1_642_3");

    write("GIDS", &gids[0]);
    write("LIDS", &lids[0]);
    write("DOMAIN", &domain[0]);
    write("FIELD", &data[0]);

    Field gdata(5, node);
    simple_write("OneD_Data", 1, 5, &gdata[0]);

    close();

    // test reading the file
    if (node == 0)
    {
        hid_t rfile = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

        std::vector<int> d(144, 0);

        hsize_t ndims[5]  = {0};
        size_t      bytes = 0;
        H5T_class_t dt_class;
        H5LTget_dataset_info(rfile, "GIDS", ndims, &dt_class, &bytes);

        EXPECT_EQ(3, ndims[0]);
        EXPECT_EQ(2, ndims[1]);
        EXPECT_EQ(4, ndims[2]);
        EXPECT_EQ(6, ndims[3]);
        EXPECT_EQ(1, ndims[4]);

        H5LTread_dataset_int(rfile, "GIDS", &d[0]);

        int c = 0;
        for (int g = 0; g < 3; ++g)
        {
            for (int k = 0; k < 2; ++k)
            {
                for (int j = 0; j < 4; ++j)
                {
                    for (int i = 0; i < 6; ++i)
                    {
                        for (int m = 0; m < 1; ++m)
                        {
                            int gid = m + 1 * (i + 6 * (j + 4 * (k + 2 * (g))));
                            EXPECT_EQ(gid, d[c]);
                            ++c;
                        }
                    }
                }
            }
        }
        EXPECT_EQ(144, c);

        H5LTread_dataset_int(rfile, "LIDS", &d[0]);
        H5Fclose(rfile);
    }
}

#endif

//---------------------------------------------------------------------------//

TEST_F(Block_Data_Simple, simple_6x4x2)
{
    if (node != 0)
        return;

    build(6, 4, 2, 1, 1, 1);

    hid_t file_id;
    file_id = H5Fcreate("simple_642.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
                        H5P_DEFAULT);

    H5LTmake_dataset_int(file_id, "GIDS", 3, dims, &gids[0]);
    H5LTmake_dataset_int(file_id, "FIELD", 3, dims, &data[0]);

    H5Fclose(file_id);
}

//---------------------------------------------------------------------------//

TEST_F(Block_Data, 1_6x4x2_3)
{
    if (node != 0)
        return;

    build(6, 4, 2, 3, 1, 1, 1, 1);

    hid_t file_id;
    file_id = H5Fcreate("simple_1_642_3.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
                        H5P_DEFAULT);

    H5LTmake_dataset_int(file_id, "GIDS", 5, dims, &gids[0]);
    H5LTmake_dataset_int(file_id, "FIELD", 5, dims, &data[0]);

    H5Fclose(file_id);

    // now read the file
    {
        hid_t rf_id = H5Fopen("simple_1_642_3.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

        // select the dataset
        hid_t dset_id = H5Dopen(rf_id, "GIDS", H5P_DEFAULT);

        // get the dataspace id
        hid_t dspace_id = H5Dget_space(dset_id);

        // select part of the slab
        hsize_t offset[5] = {0, 1, 2, 1, 0};
        hsize_t count[5]  = {1, 1, 1, 1, 1};
        hsize_t stride[5] = {1, 1, 1, 1, 1};
        hsize_t block[5]  = {3, 1, 1, 2, 1};

        // select desired dataset
        herr_t status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET,
                                            offset, stride, count, block);

        // define memory for new dataset
        hid_t memspace = H5Screate_simple(5, block, NULL);

        // select the hyperslab for this output
        hsize_t offset_out[5] = {0, 0, 0, 0, 0};
        status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                                     offset_out, stride, count, block);


        Field data(6, 0);

        // read the field
        status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, dspace_id,
                         H5P_DEFAULT, &data[0]);
        EXPECT_GE(status, 0);

        H5Sclose(dspace_id);
        H5Sclose(memspace);
        H5Dclose(dset_id);

        H5Fclose(rf_id);

        // check the data
        int ctr = 0;
        for (int g = 0; g < 3; ++g)
        {
            for (int k = 0; k < 1; ++k)
            {
                for (int j = 0; j < 1; ++j)
                {
                    for (int i = 0; i < 2; ++i)
                    {
                        for (int m = 0; m < 1; ++m)
                        {
                            int gid = m + 1 * (
                                (i+1) + 6 * ((j+2) + 4 * (
                                                 (k+1) + 2 * (g))));
                            EXPECT_EQ(gid, data[ctr]);
                            ++ctr;

                        }
                    }
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
//                 end of tstHDF5.cc
//---------------------------------------------------------------------------//
