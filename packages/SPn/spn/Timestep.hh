//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Timestep.hh
 * \author Thomas M. Evans
 * \date   Thu Apr 03 11:23:03 2014
 * \brief  Timestep class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Timestep_hh
#define spn_Timestep_hh

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "utils/Constants.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Timestep
 * \brief Timestep container for time-dependent problems.
 */
/*!
 * \example spn/test/tstTimestep.cc
 *
 * Test of Timestep.
 */
//===========================================================================//

class Timestep
{
  private:
    // >>> DATA

    // Cycle.
    int d_cycle;

    // Current timestep.
    double d_dt;

    // Stored timesteps.
    def::Vec_Dbl d_timesteps;

  public:
    // Constructor.
    Timestep()
        : d_cycle(0)
        , d_dt(constants::huge)
    {
        /*...*/
    }

    //! Set the current timestep.
    void set(double dt)
    {
        REQUIRE(dt > 0.0);

        ++d_cycle;
        d_dt = dt;

        // store the timesteps
        d_timesteps.push_back(d_dt);
    }

    //! Get the current timestep.
    double dt() const { return d_dt; }

    //! Get the current cycle.
    int cycle() const { return d_cycle; }

    //! Get the vector of timesteps.
    const def::Vec_Dbl& timesteps() const { return d_timesteps; }
};

} // end namespace profugus

#endif // spn_Timestep_hh

//---------------------------------------------------------------------------//
//                 end of Timestep.hh
//---------------------------------------------------------------------------//
