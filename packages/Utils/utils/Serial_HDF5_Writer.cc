//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Serial_HDF5_Writer.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 24 09:48:48 2014
 * \brief  Serial_HDF5_Writer member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Serial_HDF5_Writer.hh"

namespace denovo
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Serial_HDF5_Writer::Serial_HDF5_Writer()
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
void Serial_HDF5_Writer::open(const std_string &filename,
                              File_Mode         mode,
                              int               node)
{
    Require (node < b_nodes);
    Require (mode < END_FILE_MODE);

    // Set file name, mode, and node
    b_filename = filename;
    b_mode     = mode;
    b_master   = node;

    // only open the file on node 0
    if (b_node != b_master)
        return;

    // filename c-string handle
    const char* const filenm = b_filename.c_str();

    // Open the hdf5 file; create the file if not appending
    switch (mode)
    {
        case CLOBBER:
            b_file = H5Fcreate(filenm, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            break;
        case APPEND:
            b_file = H5Fopen(filenm, H5F_ACC_RDWR, H5P_DEFAULT);
            break;
        default:
            throw profugus::assertion("File mode must be CLOBBER or APPEND.");
    }

    // Set location
    HDF5_Loc loc = {b_file, ""};
    b_loc_stack.push_back(loc);

    Ensure (b_file);
    Ensure (b_loc_stack.size() == 1);
    Ensure (b_mode != END_FILE_MODE);
}

//---------------------------------------------------------------------------//
// WRITE INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Write a double to the output file.
 *
 * \param name name of the scalar data
 * \param value data to write
 */
void Serial_HDF5_Writer::write(const std_string &name,
                               double            value)
{
    Require(!is_readonly());

    if (b_node != b_master)
        return;

    // scalar dimensions
    hsize_t dims[1] = {1};

    // write the data
    H5LTmake_dataset_double(current_loc(), name.c_str(), 1, dims, &value);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write an int to the output file.
 *
 * \param name name of the scalar data
 * \param value data to write
 */
void Serial_HDF5_Writer::write(const std_string &name,
                               int               value)
{
    Require(!is_readonly());

    if (b_node != b_master)
        return;

    // scalar dimensions
    hsize_t dims[1] = {1};

    // write the data
    H5LTmake_dataset_int(current_loc(), name.c_str(), 1, dims, &value);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write a string to the output file.
 *
 * \param name name of the scalar data
 * \param value data to write
 */
void Serial_HDF5_Writer::write(const std_string &name,
                               const std_string &value)
{
    Require(!is_readonly());

    if (b_node != b_master)
        return;

    // write the data
    H5LTmake_dataset_string(current_loc(), name.c_str(), value.c_str());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write a char to the output file.
 *
 * \param name name of the scalar data
 * \param value data to write
 */
void Serial_HDF5_Writer::write(const std_string &name,
                               char              value)
{
    Require(!is_readonly());

    if (b_node != b_master)
        return;

    // scalar dimensions
    hsize_t dims[1] = {1};

    // write the data
    H5LTmake_dataset_char(current_loc(), name.c_str(), 1, dims, &value);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write a field of doubles to the output file.
 *
 * \param name name of the vector data
 * \param value data to write
 * \param size number of elements
 */
void Serial_HDF5_Writer::write(const std_string &name,
                               const double     *value,
                               std::size_t       size)
{
    Require(!is_readonly());

    if (b_node != b_master)
        return;

    // scalar dimensions
    hsize_t dims[1] = {size};

    // write the data
    H5LTmake_dataset_double(current_loc(), name.c_str(), 1, dims, value);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write a vector<int> to the output file.
 *
 * \param name name of the vector data
 * \param value data to write
 */
void Serial_HDF5_Writer::write(const std_string &name,
                               const Vec_Int    &value)
{
    Require(!is_readonly());

    if (b_node != b_master)
        return;

    // scalar dimensions
    hsize_t dims[1] = {value.size()};

    // write the data
    H5LTmake_dataset_int(current_loc(), name.c_str(), 1, dims, &value[0]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write a field of chars to the output file.
 *
 * \param name name of the vector data
 * \param value data to write
 * \param size number of elements
 */
void Serial_HDF5_Writer::write(const std_string &name,
                               const char       *value,
                               std::size_t       size)
{
    Require(!is_readonly());

    if (b_node != b_master)
        return;

    // scalar dimensions
    hsize_t dims[1] = {size};

    // write the data
    H5LTmake_dataset_char(current_loc(), name.c_str(), 1, dims, value);
}

} // end namespace denovo

//---------------------------------------------------------------------------//
//                 end of Serial_HDF5_Writer.cc
//---------------------------------------------------------------------------//
