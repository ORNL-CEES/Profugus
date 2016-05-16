//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/comm/Timing.hh
 * \author Thomas M. Evans
 * \date   Fri Nov  9 09:41:55 2007
 * \brief  Timing class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_comm_Timing_hh
#define Utils_comm_Timing_hh

#include <Utils/config.h>
#include "Timing_Diagnostics.hh"

//---------------------------------------------------------------------------//
/*!
 * \page comm_timing Macros for timing
 *
 * Four macros are defined here; these macros insert timer calls into code if
 * a global definition, UTILS_TIMING, is greater then 0.  The default value
 * is to set UTILS_TIMING == 1. They use the profugus::Timer and
 * profugus::Timing_Diagnostics classes.
 *
 * The build system sets UTILS_TIMING through the configure option \c
 * --with-timing-diagnostics.  The following settings apply:
 * - 0 turns off all TIMER macros
 * - 1 turns on TIMER, TIMER_START, TIMER_STOP
 * - 2 turns on TIMER_2, TIMER_START_2, TIMER_STOP_2 in addition to level 1
 * - 3 turns on TIMER_3, TIMER_START_3, TIMER_STOP_3 in addition to levels 1,2
 * .
 * The default is 1.
 *
 * In code,
 * \code
 * #include "comm/Timing.hh"
 *
 * TIMER( foo);
 * TIMER_START( foo);
 * // ...
 * // code interval to time
 * // ...
 * TIMER_STOP( foo);
 * TIMER_RECORD( "Snippet", foo);
 * \endcode
 * The key "Snippet" can be used to access the stored time through the
 * Timer_Diagnostics class:
 * \code
 * #ifdef UTILS_TIMING_ON
 *   vector<string> keys = Timing_Diagnostics::timer_keys();
 *   for (int i = 0; i < keys.size(); i++)
 *   {
 *       cout << keys[i] << "\t" << Timing_Diagnostics::timer_value(keys[i])
 *            << endl;
 *   }
     Timing_Diagnostics::reset_timers();
 * #endif
 * \endcode
 */

/*!
 * \def TIMER( timer_name)
 *
 * If UTILS_TIMING_ON is defined, TIMER( timer_name) expands to:
 * \code
 *     profugus::Timer timer_name
 * \endcode
 * Otherwise it is empty.
 */

/*!
 * \def TIMER_START( timer_name)
 *
 * If UTILS_TIMING > 0 TIMER_START( timer_name) expands to:
 * \code
 *     timer_name.start()
 * \endcode
 * Otherwise it is empty.
 */

/*!
 * \def TIMER_STOP( timer_name)
 *
 * If UTILS_TIMING_ON > 0, TIMER_STOP( timer_name) expands to:
 * \code
 *     timer_name.stop()
 * \endcode
 * Otherwise it is empty.
 */

/*!
 * \def TIMER_RECORD( name, timer)
 *
 * If UTILS_TIMING_ON > 0, TIMER_RECORD( name, timer) expands to:
 * \code
 *     profugus::Timing_Diagnostics::update_timer(name, timer.wall_clock())
 * \endcode
 * Otherwise it is empty.
 */
//---------------------------------------------------------------------------//

#if !defined(UTILS_TIMING)
#define UTILS_TIMING 1
#endif

//---------------------------------------------------------------------------//
/*
 * All timing operations are inactive.
 */
#if UTILS_TIMING == 0

#define TIMER( timer)

#define TIMER_START( timer)

#define TIMER_STOP( timer)

#define TIMER_RECORD( name, timer)

#define SCOPED_TIMER(name)

#endif

//---------------------------------------------------------------------------//
/*
 * Turn on basic timing operations depending on level.
 *
 * TIMER_CLOCK (defined in Timer.hh) selects "wall_clock" when MPI is on, or
 * "user_cpu" when MPI is off.
 */

// LEVEL 1
#if UTILS_TIMING > 0

#include "Timer.hh"
#include "Scoped_Timer.hh"

#define UTILS_TIMING_ON

#define TIMER( timer) profugus::Timer timer

#define TIMER_START( timer) timer.start()

#define TIMER_STOP( timer) timer.stop()

#define TIMER_RECORD( name, timer)                                      \
    profugus::Timing_Diagnostics::update_timer(name, timer.TIMER_CLOCK())

#define SCOPED_TIMER(name) \
    profugus::Scoped_Timer scoped_timer_(name)

#endif

// LEVEL 2
#if UTILS_TIMING > 1

#define UTILS_TIMING_2_ON

#define TIMER_2( timer) profugus::Timer timer

#define TIMER_START_2( timer) timer.start()

#define TIMER_STOP_2( timer) timer.stop()

#define TIMER_RECORD_2( name, timer)                                      \
    profugus::Timing_Diagnostics::update_timer(name, timer.TIMER_CLOCK())

#define SCOPED_TIMER_2(name) \
    profugus::Scoped_Timer scoped_timer_(name)

#else

#define TIMER_2( timer)

#define TIMER_START_2( timer)

#define TIMER_STOP_2( timer)

#define TIMER_RECORD_2( name, timer)

#define SCOPED_TIMER_2(name)

#endif

// LEVEL 3
#if UTILS_TIMING > 2

#define UTILS_TIMING_3_ON

#define TIMER_3( timer) profugus::Timer timer

#define TIMER_START_3( timer) timer.start()

#define TIMER_STOP_3( timer) timer.stop()

#define TIMER_RECORD_3( name, timer)                                      \
    profugus::Timing_Diagnostics::update_timer(name, timer.TIMER_CLOCK())

#define SCOPED_TIMER_3(name) \
    profugus::Scoped_Timer scoped_timer_(name)

#else

#define TIMER_3( timer)

#define TIMER_START_3( timer)

#define TIMER_STOP_3( timer)

#define TIMER_RECORD_3( name, timer)

#define SCOPED_TIMER_3(name)

#endif

#endif // Utils_comm_Timing_hh

//---------------------------------------------------------------------------//
//              end of comm/Timing.hh
//---------------------------------------------------------------------------//
