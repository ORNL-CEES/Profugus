//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/DBC.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 11:31:52 2008
 * \brief  DBC member definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <sstream>
#include "DBC.hh"

namespace profugus
{

//===========================================================================//
// ASSERTION CLASS MEMBERS
//===========================================================================//
/*!
 * \brief Default constructor for ds++/assertion class.
 *
 * This constructor creates a ds++ exception object.  This object is
 * derived form std::runtime_error and has identical functionality.  The
 * principal purpose for this class is to provide an exception class that
 * is specialized for nemesis.  See the notes on the overall class for more
 * details.
 *
 * \param msg The error message saved with the exception.
 */
assertion::assertion( std::string const & msg )
    :  Base( msg )
{
    /* empty */
}

//===========================================================================//
/*!
 * \brief Specialized constructor for profugus::assertion class.
 *
 * This constructor creates a ds++ exception object.  This object is
 * derived form std::runtime_error and has identical functionality.  This
 * constructor is specialized for use by nemesis DbC commands (Require,
 * Ensure, Check, and Insist).  It forms the error message from the test
 * condition and the file and line number of the DbC command.
 *
 * \param cond The expression that failed a DbC test.
 * \param file The source code file name that contains the DbC test.
 * \param line The source code line number that contains the DbC test.
 */
assertion::assertion(
        std::string const & cond,
        std::string const & file,
        int const line )
    : Base( build_message( cond, file, line ) )
{
    /* empty */
}

//===========================================================================//
/*! \brief Destructor for ds++/assertion class.
 *
 * We do not allow the destructor to throw!
 */
assertion::~assertion() throw()
{
    /* * */
}

//===========================================================================//
/*!
 * Build the error string (private member function).
 * \param cond Condition (test) that failed.
 * \param file The name of the file where the assertion was tested.
 * \param line The line number in the file where the assertion was tested.
 * \retval myMessage A string that contains the failed condition, the file and
 * the line number of the error.
 */
std::string assertion::build_message( std::string const & cond,
                                      std::string const & file,
                                      int         const line ) const
{
    std::ostringstream myMessage;
    myMessage << "Assertion: "
              << cond
              << ", failed in "
              << file
              << ":"
              << line;
    return myMessage.str();
}

//===========================================================================//
// NOT_IMPLEMENTED_ERROR CLASS MEMBERS
//===========================================================================//
not_implemented_error::not_implemented_error(
        const char* const msg,
        const char* const file,
        const int         line )
  : Base(build_notimpl_message(msg, file, line))
{
    /* empty */
}

//===========================================================================//
std::string not_implemented_error::build_notimpl_message(
        const char* const msg,
        const char* const file,
        const int         line)
{
    std::ostringstream out;
    out << "Regrettably, " << msg << " is not currently implemented. ";
#if UTILS_DBC > 0
    out << "\n ^^^ at " << file << ":" << line << "\n";
#else
    // Debug mode is off, so don't trouble the user with the particulars
    (void)sizeof(file);
    (void)sizeof(line);
#endif
    out << "Aborting...";

    return out.str();
}

//===========================================================================//
// validation_error CLASS MEMBERS
//===========================================================================//
validation_error::validation_error(
        const std::string& msg)
  : Base(msg)
{
    /* empty */
}


//===========================================================================//
// FREE FUNCTIONS
//===========================================================================//
/*!
 * \brief Throw a profugus::assertion for Require, Check, Ensure macros.
 * \return Throws an assertion.
 * \note We do not provide unit tests for functions whose purpose is to throw
 * or exit.
 */
void
toss_cookies(char const * const cond,
             char const * const file,
             int  const line )
{
    throw assertion( cond, file, line );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Throw a profugus::assertion for Validate macros.
 * \return Throws an assertion.
 * \note We do not provide unit tests for functions whose purpose is to throw
 * or exit.
 */
void toss_validation_cookies(
        const std::string& msg,
        const char*        file,
        int                line)
{
    std::ostringstream out;
    out << msg;
#if UTILS_DBC > 0
    out << "\n ^^^ at " << file << ":" << line << "\n";
#else
    // Debug mode is off, so don't trouble the user with the particulars
    (void)sizeof(file);
    (void)sizeof(line);
#endif

    throw validation_error(out.str());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Throw a profugus::assertion for Insist macros.
 */
void insist( const char* const  cond,
             const std::string& msg,
             const char* const  file,
             const int          line)
{
    std::ostringstream myMessage;
    myMessage <<  "Insist: " << cond << ", failed in "
              << file << ":" << line << std::endl
              << "The following message was provided:" << std::endl
              << "\"" << msg << "\"" << std::endl;
    throw assertion( myMessage.str() );
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of DBC.cc
//---------------------------------------------------------------------------//
