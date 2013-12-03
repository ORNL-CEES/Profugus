//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/DBC.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 11:31:52 2008
 * \brief  DBC class and type definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef harness_DBC_hh
#define harness_DBC_hh

#include <stdexcept>
#include <string>
#include <sstream>

// add the harness package configure
#include <Utils/config.h>

namespace nemesis
{

//===========================================================================//
/*!
 * \class assertion
 * \brief Exception notification class for Nemesis specific assertions.
 *
 * This class is derived from std::runtime_error.  In fact, this class
 * provides no significant change in functionality from std::runtime_error.
 * This class provides the following features in addition to those found in
 * std::runtime_error:
 *
 * -# nemesis::assertion does provide an alternate constructor that allows us
 *    to automatically insert file name and line location into error messages.
 * -# It provides a specialized form of std::runtime_error.  This allows
 *    nemesis code to handle nemesis specific assertions differently from
 *    generic C++ or STL exceptions.  For example
 *    \code
 *    try
 *    {
 *       throw nemesis::assertion( "My error message." );
 *    }
 *    catch ( nemesis::assertion &a )
 *    {
 *       // Catch nemesis exceptions first.
 *       cout << a.what() << endl;
 *       exit(1);
 *    }
 *    catch ( std::runtime_error &e )
 *    {
 *       // Catch general runtime_errors next
 *       cout << e.what() << endl;
 *    }
 *    catch ( ... )
 *    {
 *       // Catch everything else
 *        exit(1);
 *    }
 *    \endcode
 *
 * \note Assertion should always be thrown as objects on the stack and caught
 *       as references.
 *
 * \sa \ref Nemesis_DBC
 */
/*!
 * \example harness/test/tstDBC.cc
 *
 * Assertion and DBC examples.
 */
//===========================================================================//

class assertion : public std::logic_error
{
    typedef std::logic_error Base;
  public:
    // Default constructor for ds++/assertion class.
    explicit assertion( std::string const & msg );

    // \brief Specialized constructor for nemesis::assertion class.
    assertion( std::string const & cond,
               std::string const & file,
               int const line );

    /*! \brief Destructor for ds++/assertion class.
     * We do not allow the destructor to throw! */
    virtual ~assertion() throw();

  private:
    //! Helper function to build error message
    std::string build_message( std::string const & cond,
                               std::string const & file,
                               int         const   line ) const;
};

//---------------------------------------------------------------------------//
/*!
 * \class not_implemented_error
 * \brief Specialization for un-implemented features.
 *
 * Providing a thin specialization means this assertion can be handled
 * separately, perhaps with an extra apology to the user.
 */
class not_implemented_error : public assertion
{
    typedef assertion Base;
  public:
    not_implemented_error(
            const char* const msg,
            const char* const file,
            const int         line );

  private:
    //! Build the error message
    std::string build_notimpl_message(
            const char* const msg,
            const char* const file,
            const int         line);
};

//---------------------------------------------------------------------------//
/*!
 * \class validation_error
 * \brief Specialization for value validation
 *
 * Providing a thin specialization means this assertion can be handled
 * separately (e.g. raising a ValueError in the Python wrapper)
 */
class validation_error : public assertion
{
    typedef assertion Base;
  public:
    explicit validation_error(const std::string& msg);
};

//---------------------------------------------------------------------------//
// FREE NAMESPACE FUNCTIONS
//---------------------------------------------------------------------------//

//! Throw a nemesis::assertion for Require, Check, Ensure.
void toss_cookies(char const * const cond,
                  char const * const file,
                  int          const line );

//! Throw a nemesis::assertion for Insist.
void insist( const char* const  cond,
             const std::string& msg,
             const char* const  file,
             const int          line);

//! Raise a validation failure message
void toss_validation_cookies(
        const std::string& msg,
        const char*        file,
        int                line);

} // end namespace nemesis

//---------------------------------------------------------------------------//
// Load macros
#include "DBC_def.hh"
//---------------------------------------------------------------------------//

#endif // harness_DBC_hh

//---------------------------------------------------------------------------//
//              end of harness/DBC.hh
//---------------------------------------------------------------------------//
