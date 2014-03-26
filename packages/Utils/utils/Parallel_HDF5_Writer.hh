//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Parallel_HDF5_Writer.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 24 09:48:40 2014
 * \brief  Parallel_HDF5_Writer class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Parallel_HDF5_Writer_hh
#define utils_Parallel_HDF5_Writer_hh

#include "HDF5_IO.hh"

#ifdef H5_HAVE_PARALLEL

namespace profugus
{

//===========================================================================//
/*!
 * \class Parallel_HDF5_Writer
 * \brief Write hdf5 output in parallel (collective IO).
 */
/*!
 * \example io_utils/test/tstParallel_HDF5_Writer.cc
 *
 * Test of Parallel_HDF5_Writer.
 */
//===========================================================================//

class Parallel_HDF5_Writer : public HDF5_IO
{
    // Base class typedef.
    typedef HDF5_IO Base;

  public:
    // Constructor.
    Parallel_HDF5_Writer();

    // Open file for writing.
    void open(const std_string &filename, File_Mode mode = CLOBBER);

    // Write an integer field in parallel.
    void write(const std_string &name, const Decomp &d, const int *field);

    // Write a double field in parallel.
    void write(const std_string &name, const Decomp &d, const double *field);

  private:
    // >>> IMPLEMENTATION

    // Write the field.
    template<class T>
    void write_impl(const std_string &name, const Decomp &d, const T *field,
                    hid_t type);
};

} // end namespace profugus

#endif // H5_HAVE_PARALLEL
#endif // utils_Parallel_HDF5_Writer_hh

//---------------------------------------------------------------------------//
//                 end of Parallel_HDF5_Writer.hh
//---------------------------------------------------------------------------//
