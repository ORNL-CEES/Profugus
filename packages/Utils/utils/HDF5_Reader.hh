//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/HDF5_Reader.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 24 09:48:55 2014
 * \brief  HDF5_Reader class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_HDF5_Reader_hh
#define utils_HDF5_Reader_hh

#include "HDF5_IO.hh"

#ifdef USE_HDF5

namespace profugus
{

//===========================================================================//
/*!
 * \class HDF5_Reader
 * \brief HDF5 reader (independent IO).
 */
/*!
 * \example io_utils/test/tstHDF5_Reader.cc
 *
 * Test of HDF5_Reader.
 */
//===========================================================================//

class HDF5_Reader : public HDF5_IO_Node
{
    // Base class typedef.
    typedef HDF5_IO_Node Base;

  public:
    // Constructor.
    HDF5_Reader();

    // Open file for writing.
    void open(const std_string &filename, int node = 0);

    // Get the decomposition of a dataset.
    void get_decomposition(const std_string &name, Decomp &d);

    // Read a chunk of data from an integer data set.
    void read(const std_string &name, const Decomp &d, int *data);

    // Read a chunk of data from a double data set.
    void read(const std_string &name, const Decomp &d, double *data);

    // >>> SIMPLE READ INTERFACE

    // Read a double scalar from the output file.
    void read(const std_string &name, double &value);

    // Read an int scalar from the output file.
    void read(const std_string &name, int &value);

    // Read a string from the output file.
    void read(const std_string &name, std_string &value);

    // Read a char scalar from the output file.
    void read(const std_string &name, char &value);

    // Read a vector of doubles from the output file.
    void read(const std_string &name, Vec_Dbl &value);

    // Read a vector of ints from the output file.
    void read(const std_string &name, Vec_Int &value);

    // Read a vector of chars from the output file.
    void read(const std_string &name, Vec_Char &value);

  private:
    // >>> IMPLEMENTATION

    // Read the field.
    template<class T>
    void read_impl(const std_string &name, const Decomp &d, T *field,
                   hid_t type);
};

} // end namespace profugus

#endif // USE_HDF5
#endif // utils_HDF5_Reader_hh

//---------------------------------------------------------------------------//
//                 end of HDF5_Reader.hh
//---------------------------------------------------------------------------//
