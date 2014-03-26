//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/HDF5_IO.hh
 * \author Thomas M. Evans
 * \date   Thu Nov 07 21:53:44 2013
 * \brief  HDF5_IO class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_HDF5_IO_hh
#define utils_HDF5_IO_hh

#include <Utils/config.h>
#ifdef USE_HDF5

#include <string>
#include <vector>
#include <hdf5.h>
#include <hdf5_hl.h>

#include "harness/DBC.hh"
#include "utils/Definitions.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \page hdf5_io_page HDF5 IO Services
 *
 * There are 3 basic HDF5 IO classes that can be used to build
 * application-specific HDF5 IO services.  Two classes are used to write files
 * while a single class is used to read files.  HDF5 IO can be classified
 * under two idioms, \e collective and \e independent.  Collective IO occurs
 * when every process writes to a single HDF5 instance (not necessarily at the
 * same time).  Independent IO occurs when each process writes/reads
 * independently of the other processes.  Examples of collective operations
 * are:
 * - writing a distributed field to a single HDF5 dataset
 * - writing multiple datasets to a single HDF5 file concurrently
 * .
 * Independent operational examples are:
 * - writing datasets to multiple files concurrently
 * - writing datasets on a subset of processes non-concurrently
 * .
 * The HDF5 classes that write data are profugus::Parallel_HDF5_Writer and
 * profugus::Serial_HDF5_Writer.  Parallel_HDF5_Writer is designed for
 * collective write operations in which distributed data is written to a
 * single dataset from multiple processes.  Serial_HDF5_Writer is designed for
 * any type of independent writing; although the most common use-cases are
 * a.) writing global data from a specified IO process
 * b.) writing distributed data to multiple files concurrently
 * .
 * We do not currently support the collective operations of writing global
 * data concurrently for multiple processes to a single dataset in
 * Parallel_HDF5_Writer.
 *
 * For reading the class HDF5_Reader is provided.  This class enables \e
 * independent IO reading. The expected use-case is that every process will
 * open the appropriate HDF5 file and read some portion of data.  Because
 * there are no collisions in read operations, it is safe for multiple
 * processes to read from a single file concurrently.  Thus, the notion of
 * collective and independent reads has a different context than writing.
 * This will become clear in the simple examples in
 * test/tstParallel_HDF5_Writer.cc, test/tstSerial_HDF5_Writer.cc, and
 * test/tstHDF5_Reader.cc.
 *
 * All of the read/write classes share a simple interface to add hierarchy to
 * the HDF5 data file.  The following simple example illustrates this:
 * \code
   Parallel_HDF5_Writer w;

   // the 3-dimensional block array holds the local indices of the data-block
   // in (i,j,k) parallel space; for a 5x4x2 decomposition the I has range \c
   // [0,4] the J has range [0,3] and the K has range [0,1]

   // make a decomposition object giving the decomposition of the data field
   // we intend to write
   Parallel_HDF5_Writer::Decomp d(3); // decomposed over a 3-dim mesh

   // the global dimensions of the mesh are 100x100x100
   d.global[I] = 100;
   d.global[J] = 100;
   d.global[K] = 100;

   // the local dimensions are (for a 5x4x2 decomposition)
   d.local[I] = 20;
   d.local[J] = 25;
   d.local[K] = 50;

   // the offset on each domain
   d.offset[I] = d.local[I] * block[I];
   d.offset[J] = d.local[J] * block[J];
   d.offset[K] = d.local[K] * block[K];

   // by default the data is assumed to ordered column-major, but we will
   // explicitly add it (setting d.order = ROW_MAJOR changes the default)
   d.order = Parallel_HDF5_Writer::COLUMN_MAJOR;

   // the data has size (20*25*50) and is ordered COLUMN_MAJOR
   //   [i + Ni * (j + Nj * (k))

   // open the file for writing (clobbering an existing file)
   w.open("output.h5");

   // we want to put this data in the sub-directory (group) "distributed_data"
   w.begin_group("distributed_data");

   // write the field collectively
   w.write("data_field", d, &data[0]);

   // end this group (returning to "/" in the file hierarchy)
   w.end_group();

   // close the file
   w.close();
 * \endcode
 *
 * To read this file the following code can be used:
 *
 * Both Parallel_HDF5_Writer and Serial_HDF5_Writer will produce HDF5 files
 * that are indistinguishable.  For example, writing a distributed field using
 * the parallel writer will produce identical output to that of the serial
 * writer if the field was manually reduced in memory to a single processor
 * and output from that processor.
 */
/*!
 * \class HDF5_IO
 * \brief Base class for HDF5 read/write classes.
 *
 * This class class provides a base for collective HDF5 IO operations.  HDF5
 * IO semantics define \e collective and \e independent processes.
 *
 * \arg \e collective processes are expected to be called on every domain in a
 * parallel context
 *
 * \arg \e independent processes are not called on every domain; they may or
 * may not be defined in a parallel (MPI) context
 *
 * All of the base functions in this class are potentially \e collective, that
 * is, they are executed on every calling process.  However, they are not
 * mandated as being collective, ie. there are no MPI calling semantics in the
 * functions of this class. The HDF5_IO_Node base class defines a the same
 * member services in an independent context.
 *
 * \sa HDF5_IO_Node, \ref hdf5_io_page
 */
//===========================================================================//

class HDF5_IO
{
  public:
    //@{
    //! Typedefs.
    typedef def::Vec_Dbl         Vec_Dbl;
    typedef def::Vec_Int         Vec_Int;
    typedef def::Vec_Char        Vec_Char;
    typedef std::string          std_string;
    typedef std::vector<hsize_t> Vec_Hsize;
    //@}

    //! Flags specifying how to open/modify file
    enum File_Mode
    {
        READ = 0,       //!< Read-only access
        CLOBBER,        //!< Write access, clobbering existing files
        APPEND,         //!< Write access, appending to existing files
        END_FILE_MODE
    };

    //! Flags indicating order of multi-D field data.
    enum Data_Order
    {
        COLUMN_MAJOR = 0, //!< FORTRAN-STYLE (default)
        ROW_MAJOR         //!< C-STYLE
    };

    //! Decomposition struct for parallel and chunk field input/output.
    struct Decomp
    {
        hsize_t    ndims;
        Vec_Hsize  global;
        Vec_Hsize  local;
        Vec_Hsize  offset;
        Data_Order order;

        explicit Decomp(int n = 0, Data_Order o = COLUMN_MAJOR);
    };

  protected:
    // >>> DATA

    // HDF5 file handle.
    hid_t b_file;

    // Filename.
    std_string b_filename;

    // Nodes.
    int b_node;
    int b_nodes;

    // File access mode
    File_Mode b_mode;

  public:
    // Constructor.
    HDF5_IO();

    // Virtual destructor.
    virtual ~HDF5_IO() = 0;

    // Close output.
    virtual void close();

    // Whether we're closed
    bool closed() const { return b_file == 0; }

    // Whether we're in read-only mode
    bool is_readonly() const { return b_mode == READ; }

    //! Check for existence of a key (alias to query)
    bool exists(const std_string &name) const { return query(name); }

    // Query for a key.
    virtual bool query(const std_string &name) const;

    // >>> HIERARCHY INTERFACE

    // Begin a new group (i.e., a subdirectory)
    virtual void begin_group(const std_string &name);

    // End the current group
    virtual void end_group();

    // Return the absolute path to the current group
    virtual std_string abspath() const;

  protected:
    // >>> IMPLEMENTATION

    // Current location in HDF5 directory hierarchy (defaults to root, '/')
    struct HDF5_Loc
    {
        hid_t      handle;
        std_string group_name;
    };

    // Stack of group locations
    std::vector<HDF5_Loc> b_loc_stack;

    // Access the current location (top of the stack)
    hid_t current_loc() const
    {
        Require(!b_loc_stack.empty());
        return b_loc_stack.back().handle;
    }
};

//===========================================================================//
/*!
 * \class HDF5_IO_Node
 * \brief Base class for HDF5 read/write classes that operate on a single
 * node indepently.
 */
//===========================================================================//

class HDF5_IO_Node : public HDF5_IO
{
    // Base class typedef.
    typedef HDF5_IO Base;

  protected:
    // >>> IMPLEMENTATION

    // Master node where all writes occur.
    int b_master;

  public:
    // Constructor.
    HDF5_IO_Node();

    // Destructor.
    virtual ~HDF5_IO_Node() = 0;

    // Close file after writing.
    void close();

    // Query for a key.
    bool query(const std_string &name) const;

    // >>> HIERARCHY INTERFACE

    // Begin a new group (i.e., a subdirectory)
    void begin_group(const std_string &name);

    // End the current group
    void end_group();

    // Return the absolute path to the current group
    std_string abspath() const;
};

} // end namespace profugus

#endif // USE_HDF5
#endif // utils_HDF5_IO_hh

//---------------------------------------------------------------------------//
//                 end of HDF5_IO.hh
//---------------------------------------------------------------------------//
