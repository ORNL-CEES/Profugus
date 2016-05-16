//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/Keff_Tally.t.hh
 * \author Thomas M. Evans
 * \date   Wed May 14 13:29:40 2014
 * \brief  Keff_Tally member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Keff_Tally_t_hh
#define MC_mc_Keff_Tally_t_hh

#include "Keff_Tally.hh"

#include "harness/DBC.hh"
#include "comm/global.hh"
#include "Definitions.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Kcode solver should construct this with initial keff estimate.
 */
template <class Geometry>
Keff_Tally<Geometry>::Keff_Tally(double     keff_init,
                                 SP_Physics physics)
    : Base(physics, true)
    , d_keff_cycle(keff_init)
{
    REQUIRE(physics);

    // set the tally name
    this->set_name("keff");

    // reset tally
    reset();
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate average keff over active cycles.
 *
 * This should really only be called after begin_active_cycles; but either way
 * it will return the average keff since either initialization of this tally or
 * the last call to begin_active_cycles.
 *
 * This is only meaningful after one complete cycle. If called before then,
 * we'll return an arbitrary value.
 */
template <class Geometry>
double Keff_Tally<Geometry>::mean() const
{
    if (d_cycle < 1)
        return -1.;

    double keff = d_keff_sum / d_cycle;

    ENSURE(keff >= 0.);
    return keff;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate variance of average keff estimate over active cycles.
 *
 * This is only meaningful after two complete cycles. If called before then,
 * we'll return an arbitrary value.
 */
template <class Geometry>
double Keff_Tally<Geometry>::variance() const
{
    if (d_cycle < 2)
        return d_keff_sum * d_keff_sum;

    double var = (d_keff_sum_sq - d_keff_sum * d_keff_sum / d_cycle)
                    / (d_cycle * (d_cycle - 1));

    ENSURE(var >= 0.);
    return var;
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Track particle and do tallying.
 *
 * This only uses the particle's weight. The accumulated tally is
 * \f[
   k_{\mbox{eff}} = wl\nu\sigma_{\mbox{f}}
 * \f]
 * where \f$l\f$ is the step-length.
 */
template <class Geometry>
void Keff_Tally<Geometry>::accumulate(double            step,
                            const Particle_t &p)
{
    REQUIRE(b_physics);

    d_keff_cycle += p.wt() * step * b_physics->total(physics::NU_FISSION, p);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Begin active cycles in a kcode calculation.
 *
 * This resets the accumulated keff statistics.
 */
template <class Geometry>
void Keff_Tally<Geometry>::begin_active_cycles()
{
    d_cycle       = 0;
    d_keff_sum    = 0.;
    d_keff_sum_sq = 0.;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Begin a new keff cycle
 *
 * This clears the current accumulated path lengths.
 */
template <class Geometry>
void Keff_Tally<Geometry>::begin_cycle()
{
    d_keff_cycle = 0.;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Mark the end of a keff cycle
 *
 * This performs a global sum across processors, because the provided
 * num_particles is the total number over all domains (blocks plus sets).
 *
 * We accumulate the sum and sum-of-squares of the keff so that we can calculate
 * averages and variances.
 */
template <class Geometry>
void Keff_Tally<Geometry>::end_cycle(double num_particles)
{
    REQUIRE(num_particles > 0.);

    // Keff estimate is total nu-sigma-f reaction rate / num particles
    d_keff_cycle /= num_particles;

    // Do a global sum (since num_particles is global)
    global_sum(d_keff_cycle);

    // Accumulate first and second moments of cycles, as well as counter
    ++d_cycle;
    d_keff_sum    += d_keff_cycle;
    d_keff_sum_sq += d_keff_cycle * d_keff_cycle;

    // Store keff estimate
    d_all_keff.push_back(d_keff_cycle);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Clear accumulated keff values
 *
 * This will still preserve the latest keff value so that existing fission sites
 * will be sampled correctly.
 */
template <class Geometry>
void Keff_Tally<Geometry>::reset()
{
    d_cycle       = 0;
    d_keff_sum    = 0.;
    d_keff_sum_sq = 0.;

    d_all_keff.clear();
}

} // end namespace profugus

#endif // MC_mc_Keff_Tally_t_hh

//---------------------------------------------------------------------------//
//                 end of Keff_Tally.t.hh
//---------------------------------------------------------------------------//
