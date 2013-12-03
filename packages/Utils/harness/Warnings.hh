//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/Warnings.hh
 * \author Thomas M. Evans
 * \date   Sun Feb 26 20:54:46 2012
 * \brief  Warnings class definition.
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef harness_Warnings_hh
#define harness_Warnings_hh

#include <deque>
#include <string>
#include <sstream>
#include <iostream>

#include <harness/config.h>

namespace nemesis
{

//===========================================================================//
/*!
 * \class Warnings
 * \brief Provides warnings container for lazy-processing.
 *
 * The Warnings class is implemented as a double-ended queue: it allows warnings
 * to be pushed onto the back and popped from the front.
 */
/*!
 * \example harness/test/tstWarnings.cc
 *
 * Test of Warnings.
 */
//===========================================================================//

class Warnings
{
  private:
    //! Underlying warnings torage type
    typedef std::deque<std::string> Container;
  public:
    //! Explicitly define size type for swig wrappers
    typedef unsigned int size_type;

  private:
    // >>> DATA

    // Warnings.
    Container d_warnings;

  public:
    // Constructor.
    Warnings();

    // Add a warning.
    void add(const std::string &w);

    //! Pop the oldest warning
    std::string pop();

    //! Clear all stored warnings
    void clear() { d_warnings.clear(); }

    //! Are we empty?
    bool empty() const { return d_warnings.empty(); }

    //! Number of stored warnings.
    size_type num_warnings() const { return d_warnings.size(); }
};

//---------------------------------------------------------------------------//
// EXTERNAL DEFINITION
//---------------------------------------------------------------------------//

namespace warn
{

extern Warnings warnings;

} // end of warn

} // end namespace nemesis

//---------------------------------------------------------------------------//
/*!
 * \page warnings_controls Warnings Usage in Denovo
 *
 * The warnings class is best used through the provided macros:
 * \code
   // simple string
   ADD_WARNING("A warning");

   // ostream
   ADD_WARNING("This value " << x << " is potentially out of range");

   // get the warnings in the order that they were emplaced
   while (!NEMESIS_WARNINGS.empty())
   {
       std::cout << NEMESIS_WARNINGS.pop() << std::endl;
   }
   \endcode
 *
 * To print warnings to \c std::cout rather than save them for later, set the
 * ENABLE_NEMESIS_IMMEDIATE_WARN CMake option.
 *
 * Storing warnings is on by default, although it is up to the client to output
 * them.  Disable warning storage by setting the ENABLE_NEMESIS_WARNINGS cmake
 * option.  The code does not have to change in either case (as long as the
 * macros are used).
 *
 * The warnings object can also be accessed manually, in which case the
 * registered warnings are \b always stored, ie:
 * \code
   nemesis::warn::warnings.add("a warning");

   // or, equivalently
   NEMESIS_WARNINGS.add("a warning");

   // NEMESIS_WARNINGS is always defined, so these statements are \b always
   // equivalent
   // whether NEMESIS_WARNINGS_ENABLED is defined or not
   \endcode
 * The macro flag \c NEMESIS_WARNINGS_ENABLED is defined when warnings are
 * enabled.
 * \code
   #include <harness/config.h>

   #ifdef NEMESIS_WARNINGS_ENABLED
   cout << "Here are some warnings" << endl;
   // ...
   #else
   cout << "Warnings turned off for this build" << endl;
   #endif
   \endcode
 */
//---------------------------------------------------------------------------//

#define NEMESIS_WARNINGS ::nemesis::warn::warnings

#ifdef NEMESIS_WARNINGS_IMMEDIATE
#define ADD_WARNING(b) do { \
        std::cout << "*** Warning: " << b << std::endl; \
    } while (0)
#elif defined(NEMESIS_WARNINGS_ENABLED)
#define ADD_WARNING(b) do { \
        std::ostringstream m_; m_ << b; NEMESIS_WARNINGS.add(m_.str()); \
    } while (0)
#else
#define ADD_WARNING(b) do { } while (0)
#endif

#endif // harness_Warnings_hh

//---------------------------------------------------------------------------//
//              end of harness/Warnings.hh
//---------------------------------------------------------------------------//
