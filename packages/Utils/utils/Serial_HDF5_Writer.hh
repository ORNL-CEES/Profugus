//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Serial_HDF5_Writer.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 24 09:48:48 2014
 * \brief  Serial_HDF5_Writer class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Serial_HDF5_Writer_hh
#define utils_Serial_HDF5_Writer_hh

#include "HDF5_IO.hh"

#ifdef USE_HDF5

namespace denovo
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
    void write(const std_string &name, const double *value,
               std::size_t size);

    // Write a vector of ints to the output file.
    void write(const std_string &name, const Vec_Int &value);

    // Write a vector of chars to output file.
    void write(const std_string &name, const Vec_Char &value)
    {
        this->write(name, (value.empty() ? NULL : &value[0]), value.size());
    }

    // Write a field of chars to the output file.
    void write(const std_string &name, const char *value,
               std::size_t size);
};

} // end namespace denovo

#endif // USE_HDF5
#endif // utils_Serial_HDF5_Writer_hh

//---------------------------------------------------------------------------//
//                 end of Serial_HDF5_Writer.hh
//---------------------------------------------------------------------------//
