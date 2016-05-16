//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/utils/Serial_HDF5_Writer.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 24 09:48:48 2014
 * \brief  Serial_HDF5_Writer class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_utils_Serial_HDF5_Writer_hh
#define Utils_utils_Serial_HDF5_Writer_hh

#include <algorithm>

#include "HDF5_IO.hh"

#ifdef USE_HDF5

namespace profugus
{

//===========================================================================//
/*!
 * \class Serial_HDF5_Writer
 * \brief Write hdf5 output serially (independent IO).
 */
/*!
 * \example io_utils/test/tstSerial_HDF5_Writer.cc
 *
 * Test of Serial_HDF5_Writer.
 */
//===========================================================================//

class Serial_HDF5_Writer : public HDF5_IO_Node
{
    // Base class typedef.
    typedef HDF5_IO_Node Base;

  public:
    // Constructor.
    Serial_HDF5_Writer();

    // Open file for writing.
    void open(const std_string &filename, File_Mode mode = CLOBBER,
              int node = 0);

    // >>> WRITE INTERFACE

    // Write a double scalar to the output file.
    void write(const std_string &name, double value);

    // Write an int scalar to the output file.
    void write(const std_string &name, int value);

    // Write a string to the output file.
    void write(const std_string &name, const std_string &value);

    // Write a char to the output file.
    void write(const std_string &name, char value);

    //! Write a vector of doubles to the output file.
    void write(const std_string &name, const Vec_Dbl &value)
    {
        this->write(name, (value.empty() ? NULL : &value[0]), value.size());
    }

    // Write a field of doubles to the output file.
    void write(const std_string &name, const double *value, std::size_t size);

    // Write a vector of ints to the output file.
    void write(const std_string &name, const Vec_Int &value);

    // Write a vector of chars to output file.
    void write(const std_string &name, const Vec_Char &value)
    {
        this->write(name, (value.empty() ? NULL : &value[0]), value.size());
    }

    // Write a field of chars to the output file.
    void write(const std_string &name, const char *value, std::size_t size);

    // >>> INCREMENTAL WRITE INTERFACE

    // Create a dataset that will be incrementally modified.
    template<class T>
    void create_incremental_dataspace(const std_string &name, const Decomp &d);

    // Incrementally write into an existing dataset.
    template<class T>
    void write_incremental_data(const std_string &name, const Decomp &d,
                                const T *data);
};

//---------------------------------------------------------------------------//
// TEMPLATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Create a dataset that will be incrementally modified.
 *
 * This function is used to create a dataspace that allows incremental
 * writing.  The HDF5_IO::Decomp is used to control the size of the
 * space. Ergo,
 * \code
   // make a 5x4x2 dimensioned set that is ordered COLUMN-MAJOR
   HDF5_IO::Decomp d(3, HDF5_IO::COLUMN_MAJOR);
   d.global[0] = 5;
   d.global[1] = 4;
   d.global[2] = 2;

   // setup to incrementally write integer data
   writer.create_incremental_dataspace<int>("data", d);

   // now write using incremental_write()
 * \endcode
 *
 */
template<class T>
void Serial_HDF5_Writer::create_incremental_dataspace(
    const std_string &name,
    const Decomp     &d)
{
    if (b_node != b_master)
        return;

    // hyperslab parameters
    Vec_Hsize global(d.global);

    // HDF5 always writes data out in ROW-MAJOR (C-style) format, thus if the
    // decomposition specifies the data layout as column major we have to flip
    // all of the dimensions so that the field gets written out in the correct
    // order, ie. in ROW-MAJOR the highest dimension rotates fastest, if the
    // data is COLUMN-MAJOR this value is in the lowest dimension
    if (d.order == COLUMN_MAJOR)
    {
        std::copy(d.global.rbegin(), d.global.rend(), global.begin());
    }

    // create the dataspace for the data
    hid_t filespace = H5Screate_simple(d.ndims, &global[0], NULL);

    // create the datasets
    hid_t dset_id = H5Dcreate(
        current_loc(), name.c_str(), HDF5_Traits<T>::get_type(),
        filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Sclose(filespace);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Incrementally write into an existing dataset.
 */
template<class T>
void Serial_HDF5_Writer::write_incremental_data(
    const std_string &name,
    const Decomp     &d,
    const T          *data)
{
    if (b_node != b_master)
        return;

    // hyperslab parameters
    Vec_Hsize local(d.local);
    Vec_Hsize offset(d.offset);
    Vec_Hsize count(d.ndims, 1);  // 1-block in each dimension
    Vec_Hsize stride(d.ndims, 1); // stride 1 element at a time

    // see create_incremental_dataspace() [above] for why we may need to flip
    // the indices
    if (d.order == COLUMN_MAJOR)
    {
        std::copy(d.local.rbegin(),  d.local.rend(),  local.begin());
        std::copy(d.offset.rbegin(), d.offset.rend(), offset.begin());
    }

    // open the dataset for writing
    hid_t dset_id  = H5Dopen(current_loc(), name.c_str(), H5P_DEFAULT);
    hid_t memspace = H5Screate_simple(d.ndims, &local[0], NULL);

    // select the hyperspace
    hid_t filespace = H5Dget_space(dset_id);
    herr_t status   = H5Sselect_hyperslab(
        filespace, H5S_SELECT_SET, &offset[0], &stride[0], &count[0],
        &local[0]);
    VALIDATE(status >= 0, "Failed to select hyperslab at " << current_loc());

    // write the data
    status = H5Dwrite(dset_id, HDF5_Traits<T>::get_type(), memspace, filespace,
                      H5P_DEFAULT, data);
    VALIDATE(status >= 0, "Failed to write data on " << b_node << " for data "
             << name << " at " << current_loc());

    // close
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
}

} // end namespace profugus

#endif // USE_HDF5
#endif // Utils_utils_Serial_HDF5_Writer_hh

//---------------------------------------------------------------------------//
//                 end of Serial_HDF5_Writer.hh
//---------------------------------------------------------------------------//
