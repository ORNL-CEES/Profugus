//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/HDF5_Reader.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 24 09:48:55 2014
 * \brief  HDF5_Reader member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <numeric>
#include <algorithm>
#include "HDF5_Reader.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// READ IMPLEMENTATION
// This is defined in the implementation (.cc) file because it is used
// internally within this class only (by the public write functions).  Thus,
// we don't need to make its definition visible to other clients.  It will be
// locally instantiated in this .o.
//---------------------------------------------------------------------------//
/*!
 * \brief Read implmentation.
 */
template<class T>
void HDF5_Reader::read_impl(const std_string &name,
                            const Decomp     &d,
                            T                *field,
                            hid_t             type)
{
    REQUIRE(field != 0);
    REQUIRE(d.ndims > 0);

    // get the dataset
    hid_t dset_id = H5Dopen(current_loc(), name.c_str(), H5P_DEFAULT);

    // get the dataspace where this set resides
    hid_t dspace_id = H5Dget_space(dset_id);

    // make the slab selection
    Vec_Hsize local(d.local);
    Vec_Hsize offset(d.offset);
    Vec_Hsize count(d.ndims, 1);  // 1-block in each dimension
    Vec_Hsize stride(d.ndims, 1); // stride 1 element at a time

    // HDF5 always writes data out in ROW-MAJOR (C-style) format, thus if the
    // decomposition specifies the data layout as column major we have to flip
    // all of the dimensions so that the field gets written out in the correct
    // order, ie. in ROW-MAJOR the highest dimension rotates fastest, if the
    // data is COLUMN-MAJOR this value is in the lowest dimension
    if (d.order == COLUMN_MAJOR)
    {
        std::copy(d.local.rbegin(),  d.local.rend(),  local.begin());
        std::copy(d.offset.rbegin(), d.offset.rend(), offset.begin());
    }
    CHECK(std::accumulate(local.begin(), local.end(), 0) <=
           std::accumulate(d.global.begin(), d.global.end(), 0));

    // select the desired data in the dataset on the file
    VALIDATE(H5Sselect_hyperslab(
                  dspace_id, H5S_SELECT_SET, &offset[0], &stride[0], &count[0],
                  &local[0]) >= 0, "Unable to select hyperslab at "
              << current_loc());

    // make memory state that data will be written into
    Vec_Hsize offset_out(d.ndims, 0);
    hid_t memspace = H5Screate_simple(d.ndims, &local[0], NULL);

    // select the hyperslab for this output
    VALIDATE(H5Sselect_hyperslab(
                  memspace, H5S_SELECT_SET, &offset_out[0], &stride[0],
                  &count[0], &local[0]) >= 0, "Unable to select hyperslab "
              << " memspace at " << current_loc());

    // read the field
    VALIDATE(H5Dread(dset_id, type, memspace, dspace_id, H5P_DEFAULT,
                      field) >= 0,
              "Unable to read data at " << current_loc());

    // close handles
    H5Sclose(dspace_id);
    H5Sclose(memspace);
    H5Dclose(dset_id);
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
HDF5_Reader::HDF5_Reader()
    : Base()
{
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Open an HDF5 file for reading.
 *
 * \param filename output file name, with extension
 */
void HDF5_Reader::open(const std_string &filename,
                       int               node)
{
    REQUIRE(node < b_nodes);

    // Set file name, mode, and node
    b_filename = filename;
    b_mode     = READ;
    b_master   = node;

    // only open the file on node 0
    if (b_node != b_master)
        return;

    // filename c-string handle
    const char* const filenm = b_filename.c_str();

    // open the file
    b_file = H5Fopen(filenm, H5F_ACC_RDONLY, H5P_DEFAULT);

    // Set location
    HDF5_Loc loc = {b_file, ""};
    b_loc_stack.push_back(loc);

    ENSURE(b_file);
    ENSURE(b_loc_stack.size() == 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the decomposition of a dataset.
 */
void HDF5_Reader::get_decomposition(const std_string &name,
                                    Decomp           &d)
{
    if (b_node != b_master)
        return;

    // clear decomposition data
    d.offset.clear();
    d.global.clear();
    d.local.clear();
    d.ndims = 0;

    // hdf5 parameters
    int         rank  = 0;
    size_t      bytes = 0;
    H5T_class_t dt_class;

    // get the rank of the data
    VALIDATE(H5LTget_dataset_ndims(current_loc(), name.c_str(), &rank) >= 0,
              "Failed to read " << name << " from " << b_filename);
    CHECK(rank > 0);

    // allocate space
    d.ndims = rank;
    d.offset.resize(d.ndims);
    d.global.resize(d.ndims);
    d.local.resize(d.ndims);

    // set everything to 0
    std::fill(d.local.begin(), d.local.end(), 0);
    std::fill(d.offset.begin(), d.offset.end(), 0);

    // read global dimensions
    Vec_Hsize dims(d.ndims);
    H5LTget_dataset_info(current_loc(), name.c_str(), &dims[0], &dt_class,
                         &bytes);

    // assign the dims to the global size assuming COLUMN-MAJOR ordering
    // (ie. we flip the data dimensions)
    std::copy(dims.rbegin(), dims.rend(), d.global.begin());

    // if the storage is row-major than assign in forward order
    if (d.order == ROW_MAJOR)
    {
        std::copy(dims.begin(), dims.end(), d.global.begin());
    }

    ENSURE(d.global.size() == d.ndims);
    ENSURE(d.local.size()  == d.ndims);
    ENSURE(d.offset.size() == d.ndims);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read a chunk of data from an integer data set.
 */
void HDF5_Reader::read(const std_string &name,
                       const Decomp     &d,
                       int              *data)
{
    if (b_node != b_master)
        return;

    // pass through to read implementation
    read_impl(name, d, data, H5T_NATIVE_INT);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read a chunk of data from an integer data set.
 */
void HDF5_Reader::read(const std_string &name,
                       const Decomp     &d,
                       double           *data)
{
    if (b_node != b_master)
        return;

    // pass through to read implementation
    read_impl(name, d, data, H5T_NATIVE_DOUBLE);
}

//---------------------------------------------------------------------------//
// SIMPLE READ INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Read a double scalar from the output file.
 */
void HDF5_Reader::read(const std_string &name,
                       double           &value)
{
    if (b_node != b_master)
        return;

    VALIDATE(H5LTread_dataset_double(current_loc(), name.c_str(), &value) >= 0,
              "Failed to read " << name << " from " << b_filename);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read an int scalar from the output file.
 */
void HDF5_Reader::read(const std_string &name,
                       int              &value)
{
    if (b_node != b_master)
        return;

      VALIDATE(H5LTread_dataset_int(current_loc(), name.c_str(), &value) >= 0,
                "Failed to read " << name << " from " << b_filename);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read a char scalar from the output file.
 */
void HDF5_Reader::read(const std_string &name,
                       char              &value)
{
    if (b_node != b_master)
        return;

      VALIDATE(H5LTread_dataset_char(current_loc(), name.c_str(), &value) >= 0,
                "Failed to read " << name << " from " << b_filename);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read a string from the output file.
 */
void HDF5_Reader::read(const std_string &name,
                       std_string       &value)
{
    if (b_node != b_master)
        return;

    // hdf5 parameters
    int         rank  = 0;
    size_t      bytes = 0;
    H5T_class_t dt_class;

    // get the ndims (should be 1 for 1-D vector)
    VALIDATE(H5LTget_dataset_ndims(current_loc(), name.c_str(), &rank) >= 0,
              "Failed to read " << name << " from " << b_filename);
    CHECK(rank == 0); // rank is 0 for strings for some reason

    // read the dimensions
    hsize_t ndims[1] = {0};
    H5LTget_dataset_info(current_loc(), name.c_str(), ndims, &dt_class, &bytes);

    // make vector of bytes to store string
    std::vector<char> b(bytes);

    // get the string
    H5LTread_dataset_string(current_loc(), name.c_str(), &b[0]);

    // set the string
    value = std_string(&b[0]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read a vector of doubles from the output file.
 */
void HDF5_Reader::read(const std_string &name,
                       Vec_Dbl          &value)
{
    if (b_node != b_master)
        return;

    // hdf5 parameters
    int         rank  = 0;
    size_t      bytes = 0;
    H5T_class_t dt_class;

    // get the ndims
    VALIDATE(H5LTget_dataset_ndims(current_loc(), name.c_str(), &rank) >= 0,
              "Failed to read " << name << " from " << b_filename);
    CHECK(rank > 0);

    // read the dimensions
    Vec_Hsize ndims(rank, 0);
    H5LTget_dataset_info(current_loc(), name.c_str(), &ndims[0], &dt_class,
                         &bytes);

    // calculate size
    int size = 1;
    for (int n = 0; n < rank; ++n)
    {
        size *= ndims[n];
    }

    // size the vector
    value.resize(size);

    // get the data
    H5LTread_dataset_double(current_loc(), name.c_str(), &value[0]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read a vector of ints from the output file.
 */
void HDF5_Reader::read(const std_string &name,
                       Vec_Int          &value)
{
    if (b_node != b_master)
        return;

    // hdf5 parameters
    int         rank  = 0;
    size_t      bytes = 0;
    H5T_class_t dt_class;

    // get the ndims (should be 1 for 1-D vector)
    VALIDATE(H5LTget_dataset_ndims(current_loc(), name.c_str(), &rank) >= 0,
              "Failed to read " << name << " from " << b_filename);
    CHECK(rank > 0);

    // read the dimensions
    Vec_Hsize ndims(rank, 0);
    H5LTget_dataset_info(current_loc(), name.c_str(), &ndims[0], &dt_class,
                         &bytes);

    // calculate size
    int size = 1;
    for (int n = 0; n < rank; ++n)
    {
        size *= ndims[n];
    }

    // size the vector
    value.resize(size);

    // get the data
    H5LTread_dataset_int(current_loc(), name.c_str(), &value[0]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read a vector of chars from the output file.
 */
void HDF5_Reader::read(const std_string &name,
                       Vec_Char          &value)
{
    if (b_node != b_master)
        return;

    // hdf5 parameters
    int         rank  = 0;
    size_t      bytes = 0;
    H5T_class_t dt_class;

    // get the ndims (should be 1 for 1-D vector)
    VALIDATE(H5LTget_dataset_ndims(current_loc(), name.c_str(), &rank) >= 0,
              "Failed to read " << name << " from " << b_filename);
    CHECK(rank > 0);

    // read the dimensions
    Vec_Hsize ndims(rank, 0);
    H5LTget_dataset_info(current_loc(), name.c_str(), &ndims[0], &dt_class,
                         &bytes);

    // calculate size
    int size = 1;
    for (int n = 0; n < rank; ++n)
    {
        size *= ndims[n];
    }

    // size the vector
    value.resize(size);

    // get the data
    H5LTread_dataset_char(current_loc(), name.c_str(), &value[0]);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of HDF5_Reader.cc
//---------------------------------------------------------------------------//
