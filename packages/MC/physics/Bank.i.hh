//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Bank.i.hh
 * \author Seth R Johnson and Thomas M. Evans
 * \date   Friday April 25 16:50:47 2014
 * \brief  Member definitions of class Bank.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Bank_i_hh
#define mc_Bank_i_hh

namespace profugus
{

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

Bank::Particle_t Bank::basic_pop()
{
    REQUIRE(!empty());
    REQUIRE(!d_particles.empty());
    REQUIRE(!d_count.empty());

    Particle_t new_particle;

    if (d_count.back() > 1)
    {
        // Make a new copy to return
        new_particle = d_particles.back();

        --d_count.back();
    }
    else
    {
        // Return our original copy
        new_particle = d_particles.back();

        d_particles.pop_back();
        d_count.pop_back();
    }

    // Also update the running total
    --d_total;

    return new_particle;
}

} // end namespace mc

#endif // mc_Bank_i_hh

//---------------------------------------------------------------------------//
//              end of mc/Bank.i.hh
//---------------------------------------------------------------------------//
