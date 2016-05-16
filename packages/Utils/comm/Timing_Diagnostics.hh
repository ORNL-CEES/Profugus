//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/comm/Timing_Diagnostics.hh
 * \author Thomas M. Evans
 * \date   Fri Nov  9 09:41:55 2007
 * \brief  Timing_Diagnostics class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_comm_Timing_Diagnostics_hh
#define Utils_comm_Timing_Diagnostics_hh

#include <string>
#include <map>
#include <vector>
#include <iosfwd>

namespace profugus
{

//===========================================================================//
/*!
 * \class Timing_Diagnostics
 * \brief Class to hold timing results for diagnostic output.
 *
 * This class provides a simple interface to store timing data during a
 * simulation.  Generally, one adds a timer label and value using the
 * update_timer() function:
 * \code
 *   // start a timer
 *   ...
 *   // do work
 *   ...
 *   // stop timer
 *   profugus::Timing_Diagnostics::update_timer("Solver", time);
 * \endcode
 * There is no need to "add" the timer entry for "Solver" before an update.
 * If the key "Solver" does not exist it is added with value 0.0 before
 * applying the value.
 *
 * The easiest way to use this class is through the TIMER macros.
 */
/*!
 * \example comm/test/tstTiming.cc
 *
 * Timing diagnostics test.
 */
//===========================================================================//

class Timing_Diagnostics
{
  public:
    // Useful typedef.
    typedef std::vector<std::string> Vec_Keys;
    typedef std::map<std::string, double> Map_Timers;

  private:
    // >>> PRIVATE DATA MEMBERS

    //! Map of timers.
    static Map_Timers timers;

  public:
    // >>> FUNCTIONAL INTERFACE

    // Add a value to the timer with name key.
    static void update_timer(const std::string &key, double value);

    //! Get a timer's value.  Adds timer with name key to map.
    static double timer_value(const std::string &k) { return timers[k]; }

    //! Get number of timers in map.
    static int num_timers() { return timers.size(); }

    // Return a vector of timer keys.
    static Vec_Keys timer_keys();

    // Access timers
    static const Map_Timers& view_timers() { return timers; }

    // Reset a timer.
    static void reset_timer(const std::string &key);

    // Reset all timers from the map of timers.
    static void reset_timers();

    // Delete a timer from the map of timers.
    static void delete_timer(const std::string &key);

    // Delete all timers from the map of timers.
    static void delete_timers();

    // Print a report of the timing results
    static void report(std::ostream& os, double total_time);

  private:
    // >>> IMPLEMENTATION

    // This class is never constructed.
    Timing_Diagnostics();

    // This class is also never destructed.
    ~Timing_Diagnostics();
};

} // end namespace profugus

#endif // Utils_comm_Timing_Diagnostics_hh

//---------------------------------------------------------------------------//
//                 end of Timing_Diagnostics.hh
//---------------------------------------------------------------------------//
