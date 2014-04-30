//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Timing_Diagnostics.cc
 * \author Thomas M. Evans
 * \date   Fri Nov  9 09:41:55 2007
 * \brief  Timing_Diagnostics member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Timing_Diagnostics.hh"

#include <iostream>
#include <iomanip>
#include <Utils/config.h>

#include "harness/DBC.hh"
#include "global.hh"

namespace profugus
{
//---------------------------------------------------------------------------//
// STATIC PUBLIC FUNCTIONAL INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Add to a specified timer.
 *
 * This functions adds value to the timer with name key.  The first time this
 * is called the value is added to zero.  The timers are all static, so
 * multiple calls to the same key will keep a running tally.  To reset, call
 * reset_timer().
 *
 * Calling this function adds the timer with name key to the map of timers.
 */
void Timing_Diagnostics::update_timer(const std::string &key,
                                      double             value)
{
    timers[key] += value;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a vector of timer keys.
 */
Timing_Diagnostics::Vec_Keys Timing_Diagnostics::timer_keys()
{
    using std::string;
    using std::map;

    // keys
    Vec_Keys keys(timers.size());

    // iterators
    Vec_Keys::iterator v            = keys.begin();
    map<string, double>::iterator m = timers.begin();

    // add keys to the vector
    for (; m != timers.end(); m++, v++)
    {
        Check (v != keys.end());
        *v = m->first;
    }

    // return the vector
    return keys;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reset a timer to zero.
 *
 * Calling this function adds the timer with name key to the map of timers.
 */
void Timing_Diagnostics::reset_timer(const std::string &key)
{
    timers[key] = 0.0;
    Ensure (timers[key] == 0.0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reset all timers in the map of timers to zero.
 */
void Timing_Diagnostics::reset_timers()
{
    // iterator to timers
    std::map<std::string, double>::iterator m = timers.begin();

    // reset each timer
    for (; m != timers.end(); m++)
        m->second = 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Removes a timer with name key from the map of timers.
 */
void Timing_Diagnostics::delete_timer(const std::string &key)
{
    timers.erase(key);
    Ensure (timers.count(key) == 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Delete all timers from the map.
 */
void Timing_Diagnostics::delete_timers()
{
    // null map
    std::map<std::string, double> null;

    // swap it with timers
    timers.swap(null);
    Ensure (timers.empty());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do timing output.
 *
 * This function does global reductions across all space to get the timing
 * max/mins across processors, so it should be called from all processors.
 * The output is only processed on the input processor.
 *
 * \param out Ostream for output
 * \param total total time that the invidual timers are normalized to
 */
void Timing_Diagnostics::report(std::ostream& out, double total)
{
#if UTILS_TIMING > 0

    using std::vector;
    using std::setw;
    using std::endl;
    using std::ios;

    out.precision(4);

    // get the keys
    Vec_Keys keys = Timing_Diagnostics::timer_keys();
    const int num_keys = keys.size();

    // determine max/min
    vector<double> max(keys.size(), 0.0);
    vector<double> min(keys.size(), 0.0);

    double v = 0.0;
    for (int i = 0; i < num_keys; ++i)
    {
        // get the value of the timer
        v = Timing_Diagnostics::timer_value(keys[i]);

        // get the local value
        max[i] = v;
        min[i] = v;
    }

    // get the max/min value across all processors
    global_min(&min[0], num_keys);
    global_max(&max[0], num_keys);

    global_barrier();
    if (node() == 0)
    {
        out << endl;

        out << "===================" << endl;
        out << "Final Timing Report" << endl;
        out << "===================" << endl << endl;

        // output the max/mins
        out << setw(60) << "Routine" << setw(15) << "Max Fraction"
            << setw(15) << "Min Fraction" << endl;
        out << "==============================================="
            << "==========================================="
            << endl;
        for (int i = 0; i < num_keys; ++i)
        {
            out << setw(60) <<  keys[i] << setw(15)
                << std::scientific << max[i] / total
                << setw(15) << min[i] / total << endl;
        }
        out << "==============================================="
            << "==========================================="
            << endl;
    }
    global_barrier();

    // reset the timers
    reset_timers();
#endif
}

//---------------------------------------------------------------------------//
// PRIVATE STATIC CLASS-MEMBER DEFINITIONS
//---------------------------------------------------------------------------//

std::map<std::string, double> Timing_Diagnostics::timers;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Timing_Diagnostics.cc
//---------------------------------------------------------------------------//
