//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/comm/Scoped_Timer.hh
 * \author Seth R Johnson
 * \date   Fri Oct 11 07:50:08 2013
 * \brief  Scoped_Timer class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_comm_Scoped_Timer_hh
#define Utils_comm_Scoped_Timer_hh

#include <string>

#include "Timing_Diagnostics.hh"
#include "Timer.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Scoped_Timer
 * \brief Start a timer when constructed, stop when destructor is called.
 *
 * This is useful for timing individual functions:
 \code
    void Sweeper::sweep()
    {
        Scoped_Timer t("Sweeper.sweep");
        // ...
        // and that's it! The timer will be accumulated at the end of the
        // function.
    }

 \endcode
 */
//===========================================================================//

class Scoped_Timer
{
  private:
    // >>> DATA

    // Name of this timer
    std::string d_name;

    // Underlying timer
    Timer       d_timer;

  public:

    //! Construct with timer name and start timing
    explicit Scoped_Timer(const std::string& name)
        : d_name(name)
    {
        // Begin timing
        d_timer.start();
    }

    // When destroyed, add time to diagnostics
    ~Scoped_Timer()
    {
        d_timer.stop();
        Timing_Diagnostics::update_timer(d_name, d_timer.TIMER_CLOCK());
    }
};

} // end namespace profugus

#endif // Utils_comm_Scoped_Timer_hh

//---------------------------------------------------------------------------//
//                 end of Scoped_Timer.hh
//---------------------------------------------------------------------------//
