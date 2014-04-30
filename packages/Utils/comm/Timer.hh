//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Timer.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:27:02 2008
 * \brief  Timer class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef comm_Timer_hh
#define comm_Timer_hh

#include <Utils/config.h>

#include <iostream>
#include "harness/DBC.hh"
#include "Functions.hh"

//---------------------------------------------------------------------------//
// CONFIGURATION_DEPENDENT TIMING MACROS

#ifdef COMM_MPI
//! With MPI, report wall clock time.
#define TIMER_CLOCK wall_clock
#else
//! Without MPI, report CPU time.
#define TIMER_CLOCK user_cpu
#endif

//---------------------------------------------------------------------------//

#if defined(HAVE_SYS_TIMES_H)
#include <sys/times.h>

namespace profugus
{

//===========================================================================//
/*!
 * \class Timer
 *
 * \brief POSIX standard timer.
 *
 * The Timer class is used to calculate wall clock, user cpu, and system cpu
 * timings.  It uses the POSIX standard times function, so it should work
 * well on all (POSIX) systems.
 *
 * The POSIX implementation of timers are described in Sec.8.15 of "Advanced
 * Programming in the UNIX Environment" by Stevens.
 *
 * Usage:
 * \code
 * #include <iostream>
 * #include "comm/Timer.hh"
 * using profugus::Timer;
 *
 * Timer t;
 * t.start();
 * // do stuff
 * t.stop();
 * std::cout << t.wall_clock() << std::endl;
 * \endcode
 *
 * \example comm/test/tstTime.cc
 */
//===========================================================================//

class Timer
{
  private:
    //! Beginning wall clock time.
    double begin;

    //! Ending wall clock time.
    double end;

    //! POSIX tms structure for beginning time.
    tms tms_begin;

    //! POSIX tms structure for ending time.
    tms tms_end;

    //! The number of clock ticks per second.
    //! \sa man times
    int const posix_clock_ticks_per_second;

    //! Flag determining if timer is currently on.
    bool timer_on;

    //! True if we can access MPI timers.
    bool const isMPIWtimeAvailable;

    //! Sum of wall clock time over all intervals.
    double sum_wall;

    //! Sum of system clock time over all intervals.
    double sum_system;

    //! Sum of system clock time over all intervals.
    double sum_user;

    //! Number of time intervals.
    int num_intervals;

    //! Determine if MPI Wtime is available.
    bool setIsMPIWtimeAvailable() const;

  public:
    // Constructors and destructors.
    Timer();
    virtual ~Timer() { /* empty */ }
    Timer( Timer const & rhs );

    // Main interface functions.
    inline void start();
    inline void stop();
    inline double wall_clock() const;
    inline double system_cpu() const;
    inline double user_cpu()   const;
    inline double posix_err()  const;

    //! Return the wall clock time in seconds, summed over all intervals.
    double sum_wall_clock() const { Require(! timer_on); return sum_wall; }

    //! Return the system cpu time in seconds, summed over all intervals.
    double sum_system_cpu() const { Require(! timer_on); return sum_system; }

    //! Return the user cpu time in seconds, summed over all intervals.
    double sum_user_cpu() const { Require(! timer_on); return sum_user; }

    //! Return the number of time intervals used in the sums.
    int intervals() const { Require(! timer_on); return num_intervals; }

    inline void reset();
    void print( std::ostream &, int p = 2 ) const;
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
//! Set the beginning of the time interval.
void Timer::start()
{
    Require(! timer_on);
    timer_on = true;
    ++num_intervals;

    // set both begin and tms_begin.
    begin = wall_clock_time( tms_begin );
}

//---------------------------------------------------------------------------//
//! Set the end of the time interval.
void Timer::stop()
{
    Require( timer_on );

    // set both end and tms_end.
    end      = wall_clock_time( tms_end );
    timer_on = false;

    sum_wall   += wall_clock();
    sum_system += system_cpu();
    sum_user   += user_cpu();
}

//---------------------------------------------------------------------------//
//! Return the wall clock time in seconds, for the last interval.
double Timer::wall_clock() const
{
    Require(! timer_on);
    return( end - begin );
}

//---------------------------------------------------------------------------//
//! Return the system cpu time in seconds, for the last interval.
double Timer::system_cpu() const
{
    Require(! timer_on);
    return( tms_end.tms_stime - tms_begin.tms_stime )
        / static_cast<double>(posix_clock_ticks_per_second);
}

//---------------------------------------------------------------------------//
//! Return the user cpu time in seconds, for the last interval.
double Timer::user_cpu() const
{
    Require(! timer_on);
    return( tms_end.tms_utime - tms_begin.tms_utime )
        / static_cast<double>(posix_clock_ticks_per_second);
}

//---------------------------------------------------------------------------//
//! The error in the posix timings
double Timer::posix_err() const
{
    return 1.0/static_cast<double>(posix_clock_ticks_per_second);
}

//---------------------------------------------------------------------------//
//! Reset the interval sums.
void Timer::reset()
{
    Require(! timer_on);

    begin         = 0.0;
    end           = 0.0;
    timer_on      = false;
    sum_wall      = 0.0;
    sum_system    = 0.0;
    sum_user      = 0.0;
    num_intervals = 0;
    return;
}

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//

inline std::ostream& operator<<( std::ostream &out, const Timer &t )
{
    t.print( out, 2 );
    return out;
}

} // end namespace profugus

#elif defined(HAVE_WINDOWS_H)
#include <windows.h>

namespace profugus
{

//===========================================================================//
/*!
 * \class Timer
 * \brief Non-POSIX timer class implementation for windows.
 *
 * This class uses the windows QueryPerformanceCounter function calls to
 * calculate time information.
 *
 */
//===========================================================================//

class Timer
{
  private:
    // >>> DATA

    // Counts at beginning.
    long long d_begin;

    // Counst at end.
    long long d_end;

    // Performance counter frequency.
    long long d_f;

  public:
    // Constructor.
    Timer()
        : d_f(0)
    {
        LARGE_INTEGER *f = reinterpret_cast<LARGE_INTEGER*>(&d_f);
        QueryPerformanceFrequency(f);
    }

    //@!
    //! Start/stop timer.
    void start()
    {
        LARGE_INTEGER *begin = reinterpret_cast<LARGE_INTEGER*>(&d_begin);
        QueryPerformanceCounter(begin);
    }
    void stop()
    {
        LARGE_INTEGER *end = reinterpret_cast<LARGE_INTEGER*>(&d_end);
        QueryPerformanceCounter(end);
    }
    //@}

    //@{
    //! Return the time (equivalent to wall-clock for all queries).
    double wall_clock() const { return time(); }
    double system_cpu() const { return time(); }
    double user_cpu()   const { return time(); }
    //@}

  private:
    // >>> IMPLEMENTATION

    //! Get the time.
    double time() const { return static_cast<double>(d_end-d_begin)/d_f; }
};

} // end namespace profugus

#else

namespace profugus
{

//===========================================================================//
/*!
 * \class Timer
 * \brief Null-timer, effectively turns off all timing.
 */
//===========================================================================//

class Timer
{
  public:
    // Constructor.
    Timer() {/*...*/}

    // Main interface functions
    void start() {/*...*/}
    void stop() {/*...*/}
    double wall_clock() const { return 0.0; }
    double system_cpu() const { return 0.0; }
    double user_cpu()   const { return 0.0; }
};

} // end namespace profugus

#endif // HAVE_SYS_TIMES_H

#endif // comm_Timer_hh

//---------------------------------------------------------------------------//
//              end of comm/Timer.hh
//---------------------------------------------------------------------------//
