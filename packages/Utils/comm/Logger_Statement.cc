//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/comm/Logger_Statement.cc
 * \author Seth R Johnson
 * \date   Sat Feb 21 00:12:15 2015
 * \brief  Logger_Statement class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Logger_Statement.hh"

#include "global.hh"

namespace profugus
{
//---------------------------------------------------------------------------//
Logger_Statement::Logger_Statement(Vec_Ostream&& streams)
    : d_sinks(streams)
{
#ifdef REQUIRE_ON
    for (auto sink : d_sinks)
    {
        REQUIRE(sink);
    }
#endif

    if (!d_sinks.empty())
    {
        // Allocate a message stream if we're actually doing output
        d_message.reset(new osstream_t());
    }
    ENSURE(!d_sinks.empty() == static_cast<bool>(d_message));
}

//---------------------------------------------------------------------------//
Logger_Statement::~Logger_Statement() noexcept
{
    if (!d_message)
        return;

    try
    {
        // Add a trailing newline
        *d_message << "\n";

        // Get the string output
        const auto& message = d_message->str();

        // Write it to all the streams
        for (auto stream_ptr : d_sinks)
        {
            *stream_ptr << message << std::flush;
        }
    }
    catch (const std::exception& e)
    {
        if (profugus::node() == 0)
        {
            std::cerr << "An error occurred writing a log message: "
                      << e.what()
                      << std::endl;
        }
    }
}

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
// end of Profugus/comm/Logger_Statement.cc
//---------------------------------------------------------------------------//
