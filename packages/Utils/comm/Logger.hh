//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Profugus/comm/Logger.hh
 * \author Seth R Johnson
 * \date   Tue Feb 03 08:39:27 2015
 * \brief  Logger class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Profugus_comm_Logger_hh
#define Profugus_comm_Logger_hh

// Workaround for Windows, which loves to break code by #defining max, min,
// ERROR, and all sorts of other names that no one would ever possibly want to
// use.
#ifdef ERROR
#undef ERROR
#endif

// Workaround for LibMesh, because it will # define DEBUG.
#ifdef DEBUG
#undef DEBUG
#endif

#include <memory>
#include <iostream>

#include "Logger_Statement.hh"

namespace profugus
{
//---------------------------------------------------------------------------//
/*!
 * \page logging Logging messages in Denovo
 *
 * DEBUG messages are fine-grained diagnostics.
 *
 * DIAGNOSTIC messages are for useful information that may produce a good bit
 * of output.
 *
 * STATUS messages should be essentially the same at a given program point, no
 * matter what the execution. (For example, "Building cross sections" or
 * "Transporting step 1".)
 *
 * INFO messages typically contain information unique to the particular
 * problem: e.g. "Loaded 123 cross sections" or "Set default foo to 123".
 *
 * WARNING messages should be about unexpected or unusual behavior.
 *
 * ERROR messages are when something went wrong but we're trying to recover
 * from it or ignore it.
 *
 * CRITICAL messages are meant for "the big goodbye": explaining the last error
 * message this process will give. The prefix "!*!*!" is actually intercepted
 * by Omnibus and used as a trigger to abort the MPI run if one processor
 * aborts.
 *
 * Example messages:
 * \code

   log(DEBUG) << "Nuclide " << n << " has " << r << " reactions.";
   log_local(DIAGNOSTIC) << "Finished transporting on node " << node << ".";
   log(STATUS) << "Building solver...";
   log(INFO) << "Built solution vector with " << num << " unknowns.";
   log(WARNING) << "Nuclide 1001 was remapped to 1000";
   log_local(ERROR) << "Geometry error (lost particle) in history " << n;
   log(CRITICAL) << "Caught exception " << e.what() << "; aborting.";

 * \endcode
 */

//===========================================================================//
//! Enumeration for logging level.
enum Log_Level
{
    DEBUG = 0,  //!< Debugging messages
    DIAGNOSTIC, //!< Diagnostics about current program execution
    STATUS,     //!< Program execution status (what stage is beginning)
    INFO,       //!< Important informational messages
    WARNING,    //!< Warnings about unusual events
    ERROR,      //!< Something went wrong, but execution continues
    CRITICAL,   //!< Something went terribly wrong; we're aborting now! Bye!
    END_LOG_LEVEL
};

//===========================================================================//
/*!
 * \class Logger
 * \brief Global parallel logging for Profugus.
 *
 * This singleton class is generally accessed via the "log" free function.
 *
 * Currently the node ID is saved whenever the logger is instantiated (first
 * called), so if the communicator is changed, the original "master" node will
 * be the only one that logs during a global call.
 *
 * The class is designed to replace: \code

    if (node() == 0)
    {
        cout << ">>> Global message" << endl;
    }

    cout << ">>> Encountered " << n << " cowboys on node " << node() << endl;
    \endcode

 * with \code

    profugus::log() << "Global message" ;
    profugus::log_local() << "Encountered " << n << " cowboys on node "
                         << node() ;
 * \endcode
 *
 * The streams can be redirected at will by using the Logger accessor methods.
 *
 * \note The logging object returned by log() will not evalute the arguments if
 * no output will be displayed.
 */
/*!
 * \example comm/test/tstLogger.cc
 *
 * Test of Logger.
 */
//===========================================================================//

class Logger
{
  public:
    typedef Logger_Statement::ostream_t ostream_t;
    typedef std::shared_ptr<ostream_t>  SP_ostream;

  private:
    // >>> DATA

    //! Local and global minimum log levels
    Log_Level d_local_level;
    Log_Level d_global_level;

  public:
    // >>> CONFIGURATION

    // Set MINIMUM verbosity level for local log calls to be logged.
    void set_local_level(Log_Level level);

    // Set MINIMUM verbosity level for global log calls to be logged.
    void set_global_level(Log_Level level);

    // Set output stream (UNSAFE EXCEPT WITH GLOBAL OSTREAM!)
    void set(const std::string &key,
             ostream_t         *stream_ptr, //!< NON-OWNERSHIP POINTER
             Log_Level          min_level);
    // Set output stream (safe)
    void set(const std::string &key,
             SP_ostream         stream_sp,
             Log_Level          min_level);

    // Remove an output stream
    void remove(const std::string &key);

    // >>> STREAMING

    // Return a stream appropriate to the level for node-zero output
    Logger_Statement global_stream(Log_Level level);

    // Return a stream appropriate to the level for local-node output
    Logger_Statement local_stream(Log_Level level);

    // >>> STATIC METHODS

    static Logger& get();

  private:
    Logger();
    Logger(const Logger&);
    Logger& operator=(const Logger&);

    // Prefixes for debug/info/etc e.g. "***"
    static const char* const d_log_prefix[];

  private:
    //! Struct for output levels
    struct Sink
    {
        std::string name;  //!< Name of output sink
        Log_Level   level; //!< Output only if message >= this level
        SP_ostream  retained_stream; //!< SP to keep pointer alive
        ostream_t*  stream_ptr; //!< Actual stream to be used

        Sink(const std::string &n, Log_Level lev)
            : name(n)
            , level(lev)
            , retained_stream()
            , stream_ptr(nullptr)
        {
            /* * */
        }
    };

    // Set the active stream
    void activate_stream(ostream_t* stream);

  private:
    //! Node ID, saved when Logger is first called.
    int d_node;

    // Instead of doing something complicated like a sorted vector on name,
    // just have one sink for screen output, one for "log file" output
    Sink d_screen_output;
    Sink d_file_output;

    // Find the sink given this name
    Sink& find(const std::string &key);

    // Build output streams based on the given level
    void build_streams(
            Log_Level level,
            Logger_Statement::Vec_Ostream& streams) const;
};

//---------------------------------------------------------------------------//
//! Access the global logger instance
inline Logger& logger()
{
    return Logger::get();
}

//---------------------------------------------------------------------------//
//! Return an ostream for global (node zero only) messages
inline Logger_Statement log(Log_Level level = INFO)
{
    return Logger::get().global_stream(level);
}

//---------------------------------------------------------------------------//
//! Return an ostream for local messages
inline Logger_Statement log_local(Log_Level level = INFO)
{
    return Logger::get().local_stream(level);
}

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
#endif // Profugus_comm_Logger_hh

//---------------------------------------------------------------------------//
// end of Profugus/comm/Logger.hh
//---------------------------------------------------------------------------//
