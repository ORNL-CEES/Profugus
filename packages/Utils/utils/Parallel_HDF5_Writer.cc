//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/Parallel_HDF5_Writer.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 24 09:48:40 2014
 * \brief  Parallel_HDF5_Writer member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <algorithm>

#include "comm/global.hh"
#include "Parallel_HDF5_Writer.hh"

#ifdef H5_HAVE_PARALLEL

namespace profugus
{

//---------------------------------------------------------------------------//
// WRITE IMPLEMENTATION
// This is defined in the implementation (.cc) file because it is used
// internally within this class only (by the public write functions).  Thus,
// we don't need to make its definition visible to other clients.  It will be
// locally instantiated in this .o.
//---------------------------------------------------------------------------//
/*!
 * \brief Write implmentation.
 */
template<class T>
void Parallel_HDF5_Writer::write_impl(const std_string &name,
                                      const Decomp     &d,
                                      const T          *field,
                                      hid_t             type)
{
    REQUIRE(d.ndims > 0);

    // hyperslab parameters
    Vec_Hsize local(d.local);
    Vec_Hsize global(d.global);
    Vec_Hsize offset(d.offset);
    Vec_Hsize count(d.ndims, 1);  // 1-block in each dimension
    Vec_Hsize stride(d.ndims, 1); // stride 1 element at a time
    hsize_t   adims[] = {1};      // dimensions of attributes

    // HDF5 always writes data out in ROW-MAJOR (C-style) format, thus if the
    // decomposition specifies the data layout as column major we have to flip
    // all of the dimensions so that the field gets written out in the correct
    // order, ie. in ROW-MAJOR the highest dimension rotates fastest, if the
    // data is COLUMN-MAJOR this value is in the lowest dimension
    if (d.order == COLUMN_MAJOR)
    {
        std::copy(d.local.rbegin(),  d.local.rend(),  local.begin());
        std::copy(d.global.rbegin(), d.global.rend(), global.begin());
        std::copy(d.offset.rbegin(), d.offset.rend(), offset.begin());
    }

    // now, we could have a case where the domains don't divide evenly among
    // the cores; however, parallel HDF5 requires that the chunk size on each
    // domain is equal; thus, we get the maximum local domain size and we will
    // ue that to allocate a global space that is large enough to hold equal
    // sized chunks on every domain; then we write into the actual hyperslab
    // on each domain; this means that there will be regions of padded data in
    // the output file; however the HDF_Reader can be given the correct
    // hyperslab to read over

    // find minimum and maximum local sizes across all domains - we have to
    // use a vec-int here because the comm library *may* not have an
    // instantiation of hsize_t depending on the machine type
    Vec_Int min_c(local.begin(), local.end());
    Vec_Int max_c(local.begin(), local.end());
    profugus::global_max(&max_c[0], max_c.size());
    profugus::global_min(&min_c[0], min_c.size());

    // assign the local (common-collective : same on all domains) chunk size
    Vec_Hsize chunk(max_c.begin(), max_c.end());
    CHECK(chunk.size() == local.size());
    CHECK(local.size() == d.ndims);

    // make the global chunk size - this is the chunk size * number of
    // domains; first we have to determine the number of domains that can be
    // calculated by integer-dividing the true global size by the minimum
    // chunk size
    for (int n = 0; n < d.ndims; ++n)
    {
        // find the decomposition in this dimension
        int num_domains = global[n] / min_c[n];
        CHECK(global[n] <= num_domains * chunk[n]);

        // calculate the new global for equal chunk sizes
        global[n] = num_domains * chunk[n];
     }

    // create the dataspace for the Data

    // global filespace
    hid_t filespace = H5Screate_simple(d.ndims, &global[0], NULL);
    // local filespace
    hid_t memspace  = H5Screate_simple(d.ndims, &local[0], NULL);
    // attribute filespace
    hid_t attspace  = H5Screate_simple(1, adims, NULL);

    // create the datasets (this is a collective and each domain gets the same
    // chunk size)
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, d.ndims, &chunk[0]);
    hid_t dset_id = H5Dcreate(current_loc(), name.c_str(), type,
                              filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);

    // create the attribute
    hid_t aset_id = H5Acreate_by_name(
        current_loc(), name.c_str(), "data_order", H5T_NATIVE_INT, attspace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // close objects
    H5Pclose(plist_id);
    H5Sclose(filespace);
    H5Sclose(attspace);

    // select the hyperslab (this is in the "true" domain of the data)
    filespace = H5Dget_space(dset_id);
    VALIDATE(H5Sselect_hyperslab(
                  filespace, H5S_SELECT_SET, &offset[0], &stride[0],
                  &count[0], &local[0]) >= 0,
              "Failed to select hyperslab at " << current_loc());;

    // create property list for collective I/O
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // write the data
    VALIDATE(H5Dwrite(dset_id, type, memspace, filespace, plist_id,
                       &field[0]) >= 0,
              "Failed to write data on " << b_node << " at " << current_loc());

    // write the attribute
    attspace  = H5Aget_space(aset_id);
    int order = d.order;
    VALIDATE(H5Awrite(aset_id, H5T_NATIVE_INT, &order) >= 0,
              "Failed to attach attribute to " << name.c_str());

    // close
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(attspace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Aclose(aset_id);
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Parallel_HDF5_Writer::Parallel_HDF5_Writer()
    : Base()
{
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Open an HDF5 file for output.
 *
 * \param filename output file name, with extension
 * \param mode specify
 */
void Parallel_HDF5_Writer::open(const std_string &filename,
                                File_Mode         mode)
{
    REQUIRE(mode < END_FILE_MODE);

    // Set file name, mode, and node
    b_filename = filename;
    b_mode     = mode;

    // filename c-string handle
    const char* const filenm = b_filename.c_str();

    // set property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, profugus::communicator, MPI_INFO_NULL);

    // Open the hdf5 file; create the file if not appending
    switch (mode)
    {
        case CLOBBER:
            b_file = H5Fcreate(filenm, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
            break;
        case APPEND:
            b_file = H5Fopen(filenm, H5F_ACC_RDWR, plist_id);
            break;
        default:
            throw profugus::assertion("File mode must be CLOBBER or APPEND.");
    }

    // close the properties
    H5Pclose(plist_id);

    // Set location
    HDF5_Loc loc = {b_file, ""};
    b_loc_stack.push_back(loc);

    ENSURE(b_file);
    ENSURE(b_loc_stack.size() == 1);
    ENSURE(b_mode != END_FILE_MODE);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write integer field to HDF5 file in parallel.
 */
void Parallel_HDF5_Writer::write(const std_string &name,
                                 const Decomp     &d,
                                 const int        *field)
{
    // pass through to parallel write implementation
    write_impl(name, d, field, H5T_NATIVE_INT);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write double field to HDF5 file in parallel.
 */
void Parallel_HDF5_Writer::write(const std_string &name,
                                 const Decomp     &d,
                                 const double     *field)
{
    // pass through to parallel write implementation
    write_impl(name, d, field, H5T_NATIVE_DOUBLE);
}

} // end namespace profugus

#endif // H5_HAVE_PARALLEL

//---------------------------------------------------------------------------//
//                 end of Parallel_HDF5_Writer.cc
//---------------------------------------------------------------------------//
