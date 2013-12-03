//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Timer.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:27:02 2008
 * \brief  Timer member definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iomanip>
#include "Timer.hh"

namespace profugus
{

#ifdef HAVE_SYS_TIMES_H
#include <unistd.h>

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Timer::Timer()
    : begin(0.0)
    , end(0.0)
    , posix_clock_ticks_per_second( sysconf( _SC_CLK_TCK ) )
    , timer_on(false)
    , isMPIWtimeAvailable(setIsMPIWtimeAvailable())
    , sum_wall(0.0)
    , sum_system(0.0)
    , sum_user(0.0)
    , num_intervals(0)
{
    reset();
}

//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Print out a timing report.
 */
void Timer::print( std::ostream &out, int p ) const
{
    using std::setw;
    using std::endl;
    using std::ios;

    out.setf(ios::fixed, ios::floatfield);
    out.precision(p);
    out << '\n';

    if ( num_intervals > 1 )
        out << "LAST INTERVAL: " << '\n';

    out << setw(20) << "WALL CLOCK TIME: " << wall_clock() << " sec." << '\n';
    out << setw(20) << "  USER CPU TIME: " << user_cpu()   << " sec." << '\n';
    out << setw(20) << "SYSTEM CPU TIME: " << system_cpu() << " sec." << '\n';
    out << '\n';

    if ( num_intervals > 1 )
    {
        out << "OVER " << num_intervals << " INTERVALS: " << '\n';
        out << setw(20) << "WALL CLOCK TIME: " << sum_wall_clock()
            << " sec." << '\n';
        out << setw(20) << "  USER CPU TIME: " << sum_user_cpu()
            << " sec." << '\n';
        out << setw(20) << "SYSTEM CPU TIME: " << sum_system_cpu()
            << " sec." << '\n';
        out << '\n';
    }
    out.flush();
}

//---------------------------------------------------------------------------//
// Is this an MPI or Posix timer?
//---------------------------------------------------------------------------//

bool Timer::setIsMPIWtimeAvailable() const
{
#ifdef COMM_SCALAR
    return false;
#else
    return true;
#endif
}

#endif // HAVE_SYS_TIMES_H

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Timer.cc
//---------------------------------------------------------------------------//
