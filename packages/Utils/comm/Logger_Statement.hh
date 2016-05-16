//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Utils/comm/Logger_Statement.hh
 * \author Seth R Johnson
 * \date   Sat Feb 21 00:12:15 2015
 * \brief  Logger_Statement class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_comm_Logger_Statement_hh
#define Utils_comm_Logger_Statement_hh

#include <vector>
#include <sstream>
#include <iostream>
#include <memory>
#include "harness/DBC.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Logger_Statement
 * \brief Support class for Logger that emulates an ostream.
 *
 * This class is designed to intercept ostream-like output, save it to a
 * buffer, and send it to multiple streams when it reaches the end of its
 * scope.
 *
 * It should never be stored; its lifetime should be the the scope of the
 * single statement from which it's created.
 */
//===========================================================================//

class Logger_Statement
{
    typedef Logger_Statement This;
  public:
    // >>> TYPEDEFS

    //! Output string type
    typedef std::ostream ostream_t;

    //! Vector of pointers to output streams
    typedef std::vector<ostream_t*> Vec_Ostream;

    //! Function signature for a stream maniupulator (such as endl)
    typedef ostream_t& (*Stream_Manipulator)(ostream_t&);

  private:
    // >>> DATA

    //! String stream type compatible with ostream type
    typedef std::basic_ostringstream<ostream_t::char_type,
                                     ostream_t::traits_type> osstream_t;

    //! String stream for saving this log message
    std::unique_ptr<osstream_t> d_message;

    //! Vector of "sinks" to output to
    Vec_Ostream d_sinks;

  public:

    // Construct with streams to send message to
    explicit Logger_Statement(Vec_Ostream&& streams);

    // Send message on destruction
    ~Logger_Statement() noexcept;

    // Allow moving but not copying
    Logger_Statement(Logger_Statement&&) = default;
    Logger_Statement& operator=(Logger_Statement&&) = default;
    Logger_Statement(const Logger_Statement&) = delete;
    Logger_Statement& operator=(const Logger_Statement&) = delete;

    /*!
     * \brief Act like an ostream, but return ourself.
     *
     * This allows us to intelligently disable writing expensive operations to
     * the stream if they're not going to be output. If we're saving output,
     * write the given data to the string stream.
     */
    template<class T>
    This& operator<<(const T& rhs)
    {
        if (d_message)
        {
            *d_message << rhs;
        }
        return *this;
    }

    /*!
     * \brief Specialization on const char* to reduce object size.
     */
    This& operator<<(const char* rhs)
    {
        if (d_message)
        {
            *d_message << rhs;
        }
        return *this;
    }

    /*!
     * \brief Accept manipulators such as std::endl.
     *
     * This allows us to intelligently disable writing expensive operations to
     * the stream if they're not going to be output.
     */
    This& operator<<(Stream_Manipulator manip)
    {
        if (d_message)
        {
            manip(*d_message);
        }
        return *this;
    }
};

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
#endif // Utils_comm_Logger_Statement_hh

//---------------------------------------------------------------------------//
// end of Profugus/comm/Logger_Statement.hh
//---------------------------------------------------------------------------//
