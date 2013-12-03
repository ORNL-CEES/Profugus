//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/SWIG.hh
 * \author Seth R Johnson
 * \date   Fri Feb 15 10:22:01 2013
 * \brief  Swig functions and macros that get exposed to C++ and SWIG
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef harness_SWIG_hh
#define harness_SWIG_hh

/*! \def SWIG_PRIVATE
 *
 * Access control for SWIG. Use this to prevent member functions from being
 * wrapped in Python and to hide typedefs from the SWIG wrapper code. It is
 * necessary to %include this file from every SWIG module (or %import another
 * module that %includes it).
 *
 * C++ will always see this as a "public" statement, so don't use this to try to
 * hide private data from C++.
 *
 * \code
 class Foo
 {
   public:
     void do_something();

   SWIG_PRIVATE:
     // Allow construction from C++ layer but not from Python
     Foo();
 };
 \endcode
 */
#ifdef SWIG
#define SWIG_PRIVATE private
#else
#define SWIG_PRIVATE public
#endif

#endif // harness_SWIG_hh

//---------------------------------------------------------------------------//
//              end of harness/SWIG.hh
//---------------------------------------------------------------------------//
