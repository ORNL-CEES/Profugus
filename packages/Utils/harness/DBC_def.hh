//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/DBC_def.hh
 * \author Seth R Johnson
 * \date   Tue Mar 05 09:32:46 2013
 * \brief  DBC Macro definitions
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * We provide the separate Macro definition file to make dealing with SCALE
 * (which redefines our macros) easier.
 *
 * \warning NEVER INCLUDE THIS FILE DIRECTLY! It is only for use in the DBC
 *       header files in the harness.
 *
 * \note There should be no include guards in this file. For normal cases, we
 *       rely on the DB.hh guards to prevent double inclusion. Otherwise, when
 *       redefining DBC macros, we *need* to be able to include these contents
 *       twice.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \page Nemesis_DBC Using the Nemesis Design-by-Contract Macros
 *
 * \section ddbc Using the Nemesis Design-by-Contract Macros
 *
 * The assertion macros are intended to be used for validating preconditions
 * which must be true in order for following code to be correct, etc.  For
 * example,
 *
 * \code
 * Assert( x > 0. );
 * y = sqrt(x);
 * \endcode
 *
 * If the assertion fails, the code should just bomb.  Philosophically, it
 * should be used to feret out bugs in preceding code, making sure that prior
 * results are within reasonable bounds before proceeding to use those
 * results in further computation, etc.
 *
 * These macros are provided to support the Design By Contract formalism.
 * The activation of each macro is keyed off a bit in the DBC macro which can
 * be specified on the command line:
\verbatim
     Bit     DBC macro affected
     ---     ------------------
      0      Require
      1      Check
      2      Ensure
\endverbatim
 *
 * So for instance, \c -DUTILS_DBC=7 turns them all on, \c -DUTILS_DBC=0
 * turns them all off, and \c -DUTILS_DBC=1 turns on \c Require but turns
 * off \c Check and \c Ensure.  The default is to have them all enabled.
 *
 * The \c Insist macro is akin to the \c Assert macro, but it provides the
 * opportunity to specify an instructive message.  The idea here is that you
 * should use Insist for checking things which are more or less under user
 * control.  If the user makes a poor choice, we "insist" that it be
 * corrected, providing a corrective hint.
 *
 * \note We provide a way to eliminate assertions, but not insistings.  The
 * idea is that \c Assert is used to perform sanity checks during program
 * development, which you might want to eliminate during production runs for
 * performance sake.  Insist is used for things which really really must be
 * true, such as "the file must've been opened", etc.  So, use \c Assert for
 * things which you want taken out of production codes (like, the check might
 * inhibit inlining or something like that), but use Insist for those things
 * you want checked even in a production code.
 */
/*!
 * \def Require(condition)
 *
 * Pre-condition checking macro.  On when UTILS_DBC & 1 is true.
 */
/*!
 * \def Check(condition)
 *
 * Intra-scope checking macro.  On when UTILS_DBC & 2 is true.
 */
/*!
 * \def Ensure(condition)
 *
 * Post-condition checking macro.  On when UTILS_DBC & 4 is true.
 */
/*!
 * \def Remember(code)
 *
 * Add code to compilable code.  Used in the following manner:
 * \code
 *     Remember (int old = x;)
 *     // ...
 *     Ensure (x == old);
 * \endcode
 * On when UTILS_DBC & 4 is true.
 */
/*!
 * \def Insist(condition, message)
 *
 * Inviolate check macro.  Insist is always on.
 */
/*!
 * \def Not_Implemented(feature)
 *
 * Throw an error when the given feature string is not implemented. Always on.
 *
 * Note that we don't use an intermediate function to throw the error, because
 * that causes lots of "no return statement in function returning non-void"
 * errors. Better to let the compiler know that we're throwing directly.
 *
 * If UTILS_DBC is nonzero, we print out the file and line of failure.
 * Otherwise we hide it from the user.
 */
/*!
 * \def Validate(condition, message_stream)
 *
 * Throw a user-oriented verbose error when condition is not met. The
 * message_stream is passed directly into a string stream, so it can include
 * what value failed in the output. Because "condition" can be a complicated
 * piece of code, we don't echo it to the user.
 *
 * If UTILS_DBC is nonzero, we print out the file and line of failure.
 * Otherwise we hide it from the user.
 */
//---------------------------------------------------------------------------//

#if !defined(UTILS_DBC)
#define UTILS_DBC 7
#endif

/* The following definitions of assertions help ensure that they're always
 * followed by a semicolon, and that if assertions are disabled we don't get
 * "unused variable" warnings.
 *
 * See:
 * http://cnicholson.net/2009/02/stupid-c-tricks-adventures-in-assert/
 */
#define UTILS_ASSERT_(COND) \
    do { if (!(COND)) ::nemesis::toss_cookies( \
            #COND, __FILE__, __LINE__); } while (0)
#define UTILS_NOASSERT_(COND) \
    do { (void)sizeof(COND); } while (0)


#if UTILS_DBC & 1
#define REQUIRE_ON
#define Require(c) UTILS_ASSERT_(c)
#else
#define Require(c) UTILS_NOASSERT_(c)
#endif

#if UTILS_DBC & 2
#define CHECK_ON
#define Check(c) UTILS_ASSERT_(c)
#define Assert(c) UTILS_ASSERT_(c)
#else
#define Check(c) UTILS_NOASSERT_(c)
#define Assert(c) UTILS_NOASSERT_(c)
#endif

#if UTILS_DBC & 4
#define REMEMBER_ON
#define ENSURE_ON
#define Ensure(c) UTILS_ASSERT_(c)
#define Remember(c) c
#else
#define Ensure(c) UTILS_NOASSERT_(false)
#define Remember(c)
#endif

#define Insist(COND, MSG) \
    do { if (!(COND)) ::nemesis::insist( \
            #COND, MSG, __FILE__, __LINE__); } while (0)

#define Not_Implemented(MSG) \
    throw ::nemesis::not_implemented_error(MSG, __FILE__, __LINE__)

#define Validate(COND, MSG_STREAM) \
    do \
    { \
        if (!(COND)) \
        { \
            std::ostringstream msg; \
            msg << MSG_STREAM; \
            ::nemesis::toss_validation_cookies(msg.str(), __FILE__, __LINE__); \
        } \
    } while (0)


//---------------------------------------------------------------------------//
//              end of harness/DBC_def.hh
//---------------------------------------------------------------------------//
