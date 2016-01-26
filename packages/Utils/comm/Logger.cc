//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   comm/Logger.cc
 * \author Seth R Johnson
 * \date   Tue Feb 03 08:39:27 2015
 * \brief  Logger class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Logger.hh"

#include "global.hh"

namespace profugus
{
//---------------------------------------------------------------------------//
// STATIC DATA
//---------------------------------------------------------------------------//
const char* const Logger::d_log_prefix[] = {
    "",       // DEBUG
    "",       // DIAGNOSTIC
    "::: ",   // STATUS
    ">>> ",   // INFO
    "*** ",   // WARNING
    "!!! ",   // ERROR
    "!*!*! ", // CRITICAL
};

//---------------------------------------------------------------------------//
// LOGGER
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Logger::Logger()
    : d_local_level(DIAGNOSTIC)
    , d_global_level(DIAGNOSTIC)
    , d_node(profugus::node())
    , d_screen_output("screen", DEBUG)
    , d_file_output("file",     END_LOG_LEVEL)
{
    // Default screen output is cerr
    d_screen_output.stream_ptr = &std::cerr;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set verbosity level for local log calls.
 */
void Logger::set_local_level(Log_Level level)
{
    REQUIRE(level < END_LOG_LEVEL);
    d_local_level = level;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set verbosity level for global log calls.
 */
void Logger::set_global_level(Log_Level level)
{
    REQUIRE(level < END_LOG_LEVEL);
    d_global_level = level;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set an output handle
 *
 * \warning This is UNSAFE except with global ostreams such as cout!
 *
 * Since the Logger is global, it will almost certainly exceed the scope of any
 * local stream, leading to dereferencing of deallocated data.
 *
 * If you absolutely must give a raw pointer here, make *sure* to call
 * \c remove() on it when the pointer's reference is destroyed.
 */
void Logger::set(const std::string &key,
                 ostream_t         *stream_ptr,
                 Log_Level          min_level)
{
    REQUIRE(stream_ptr);
    REQUIRE(min_level < END_LOG_LEVEL);

    Sink& sink = find(key);

    sink.name = key;
    sink.level = min_level;
    sink.stream_ptr = stream_ptr;
    sink.retained_stream = SP_ostream();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set an output handle
 *
 * This is the preferred way to set a logger handle because it has
 * reference-counted semantics.
 */
void Logger::set(const std::string &key,
                 SP_ostream         stream_sp,
                 Log_Level          min_level)
{
    REQUIRE(stream_sp);
    REQUIRE(min_level < END_LOG_LEVEL);

    Sink& sink = find(key);

    sink.name = key;
    sink.level = min_level;
    sink.stream_ptr = stream_sp.get();
    sink.retained_stream = stream_sp;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Remove an output handle.
 */
void Logger::remove(const std::string &key)
{
    Sink& sink = find(key);

    sink.stream_ptr = nullptr;
    sink.retained_stream.reset();
}

//---------------------------------------------------------------------------//
// ACCESSORS
//---------------------------------------------------------------------------//
/*!
 * \brief Return a stream appropriate to the level for node-zero output
 */
Logger_Statement Logger::global_stream(Log_Level level)
{
    REQUIRE(level < END_LOG_LEVEL);

    Logger_Statement::Vec_Ostream streams;

    // Only add streams on node zero
    if (d_node == 0 && level >= d_global_level)
    {
        build_streams(level, streams);
    }

    // Create the logger statement (moving the vec streams for efficiency)
    Logger_Statement result(std::move(streams));

    // Pipe prefix to the stream before returning
    result << d_log_prefix[level];

    // Return the expiring logger_statement, implicit move
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a stream appropriate to the level for local-node output
 */
Logger_Statement Logger::local_stream(Log_Level level)
{
    REQUIRE(level < END_LOG_LEVEL);

    Logger_Statement::Vec_Ostream streams;

    if (level >= d_local_level)
    {
        build_streams(level, streams);
    }

    // Create the logger statement (moving the vec streams for efficiency)
    Logger_Statement result(std::move(streams));

    // Pipe prefix to the stream before returning
    result << d_log_prefix[level];

    // Return the expiring logger_statement, implicit move
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Remove an output handle.
 */
Logger::Sink& Logger::find(const std::string &key)
{
    if (key == "screen")
    {
        return d_screen_output;
    }
    else if (key == "file")
    {
        return d_file_output;
    }
    else
    {
        VALIDATE(false, "Currently only screen and file are supported "
                 "log keys; '" << key << "' is invalid.");
    }

    // Squelch compiler errors
    return find(key);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build output streams based on the given level
 */
void Logger::build_streams(
        Log_Level level,
        Logger_Statement::Vec_Ostream& streams) const
{
    REQUIRE(streams.empty());
    for (const Sink& s : {d_screen_output, d_file_output})
    {
        if (s.stream_ptr && (level >= s.level))
        {
            streams.push_back(s.stream_ptr);
        }
    }
}

//---------------------------------------------------------------------------//
// STATIC METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Access global logging instance.
 *
 * \warning This should be accessed after MPI is initialized because we call
 * MPI_rank.
 */
Logger& Logger::get()
{
    static Logger d_instance;
    CHECK(std::end(d_log_prefix) - std::begin(d_log_prefix)
          == END_LOG_LEVEL);
    return d_instance;
}

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
// end of Profugus/comm/Logger.cc
//---------------------------------------------------------------------------//
