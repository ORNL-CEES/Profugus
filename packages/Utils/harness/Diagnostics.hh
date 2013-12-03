//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/Diagnostics.hh
 * \author Thomas M. Evans
 * \date   Mon Feb  2 20:06:38 2009
 * \brief  Diagnostics for runtime info.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef harness_Diagnostics_hh
#define harness_Diagnostics_hh

#include <map>
#include <vector>
#include <string>

#include <Utils/config.h>

namespace profugus
{

//===========================================================================//
/*!
 * \namespace Diagnostics
 * \brief Allows code clients to register diagnostics data during runtime.
 *
 * This namespace defines maps for the following types:
 * - int
 * - double
 * - map<int, int>
 * - vector<int>
 * - vector<double>
 * .
 * The key for each map is a std::string.
 *
 * These maps can be used to store diagnostic quantities.  Because they have
 * global, file scope they can be accessed from any routine.  The general
 * usage is as follows:
 * \code
 *   nemesis::Diagnostics::integers["Num_cells_per_proc"] = 1000;
 * \endcode
 * A compile-time switching mechanism for using these maps is provided by the
 * macros DIAGNOSTICS_ONE, DIAGNOSTICS_TWO, and DIAGNOSTICS_THREE.
 */
/*!
 * \example harness/test/tstDiagnostics.cc
 *
 * description of example
 */
//===========================================================================//

namespace Diagnostics
{

//@{
//! Mapped typedefs.
typedef std::map<std::string, int>                   Diag_Integer_t;
typedef std::map<std::string, double>                Diag_Double_t;
typedef std::map< std::string, std::vector<int> >    Diag_Vec_Integer_t;
typedef std::map< std::string, std::map<int, int> >  Diag_Map_Integer_t;
typedef std::map< std::string, std::vector<double> > Diag_Vec_Double_t;
//@}

//! Map of integer data.
extern Diag_Integer_t integers;

//! Map of floating point data.
extern Diag_Double_t doubles;

//! Map of vector, integer data.
extern Diag_Vec_Integer_t vec_integers;

//! Map of mapped integer data.
extern Diag_Map_Integer_t map_integers;

//! Map of vector, double data.
extern Diag_Vec_Double_t vec_doubles;

} // end of namespace Diagnostics

} // end namespace profugus

//---------------------------------------------------------------------------//
/*!
 * \page nemesis_diagnostics Diagnostics Levels
 *
 * The diagnostics can be turned on in three different levels based on logical
 * bit comparisions.  The following shows the levels:
 * - Bit 0, (001), activates Level 1
 * - Bit 1, (010), activates Level 2
 * - Bit 2, (100), activates Level 3
 * .
 * The following integer settings activate Levels in the following way:
 * - 0 all off
 * - 1 Level 1
 * - 2 Level 2
 * - 3 Level 1, Level 2
 * - 4 Level 3
 * - 5 Level 1, Level 3
 * - 6 Level 2, Level 3
 * - 7 Level 1, Level 2, Level 3
 * .
 * Thus setting \c --with-diagnostics=7 at configure time will turn on all
 * levels.  The default setting is 1.
 *
 * The intent is to use Level 1 for high-level, low cost diagnostics that are
 * always active (ie. User "Education").  Levels 2 and 3 are for low-level
 * diagnostics that could incur a performance penalty.  However, all of these
 * usages are up to the client.
 */
/*!
 * \def DIAGNOSTICS_ONE(Diagnostics::member)
 *
 * Single-line statement macro for diagnostics Level 1:
 * \code
 *     DIAGNOSTICS_ONE(integers["Variable"] = 1);
 * \endcode
 * On when NEMESIS_DIAGNOSTICS & 1 is true.  Defines
 * NEMESIS_DIAGNOSTICS_LEVEL_1.
 */
/*!
 * \def DIAGNOSTICS_TWO(Diagnostics::member)
 *
 * Single-line statement macro for diagnostics Level 2:
 * \code
 *     DIAGNOSTICS_TWO(integers["Variable"] = 1);
 * \endcode
 * On when NEMESIS_DIAGNOSTICS & 2 is true.  Defines
 * NEMESIS_DIAGNOSTICS_LEVEL_2.
 */
/*!
 * \def DIAGNOSTICS_THREE(Diagnostics::member)
 *
 * Single-line statement macro for diagnostics Level 3:
 * \code
 *     DIAGNOSTICS_THREE(integers["Variable"] = 1);
 * \endcode
 * On when NEMESIS_DIAGNOSTICS & 4 is true.  Defines
 * NEMESIS_DIAGNOSTICS_LEVEL_3.
 */
//---------------------------------------------------------------------------//

#if !defined(NEMESIS_DIAGNOSTICS)
#define NEMESIS_DIAGNOSTICS 1
#endif

#if NEMESIS_DIAGNOSTICS & 1
#define NEMESIS_DIAGNOSTICS_LEVEL_1
#define DIAGNOSTICS_ONE(member) nemesis::Diagnostics::member
#else
#define DIAGNOSTICS_ONE(member)
#endif

#if NEMESIS_DIAGNOSTICS & 2
#define NEMESIS_DIAGNOSTICS_LEVEL_2
#define DIAGNOSTICS_TWO(member) nemesis::Diagnostics::member
#else
#define DIAGNOSTICS_TWO(member)
#endif

#if NEMESIS_DIAGNOSTICS & 4
#define NEMESIS_DIAGNOSTICS_LEVEL_3
#define DIAGNOSTICS_THREE(member) nemesis::Diagnostics::member
#else
#define DIAGNOSTICS_THREE(member)
#endif

#endif // harness_Diagnostics_hh

//---------------------------------------------------------------------------//
//              end of harness/Diagnostics.hh
//---------------------------------------------------------------------------//
