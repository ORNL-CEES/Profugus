//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Keff_Tally.hh
 * \author Thomas M. Evans
 * \date   Wed May 14 13:29:40 2014
 * \brief  Keff_Tally class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Keff_Tally_hh
#define mc_Keff_Tally_hh

#include "Tally.hh"
#include "utils/Definitions.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Keff_Tally
 * \brief Use path length to estimate eigenvalue and variance
 *
 * Tally keff during KCode operation using path-length accumulators.
 */
/*!
 * \example mc/test/tstKeff_Tally.cc
 *
 * Test of Keff_Tally.
 */
//===========================================================================//

class Keff_Tally : public Pathlength_Tally
{
    typedef Tally Base;

  public:
    //@{
    //! Misc typedefs.
    typedef def::Vec_Dbl Vec_Dbl;
    //@}

  private:
    // >>> DATA

    //! Estimate of k-effective from this cycle
    double d_keff_cycle;

    //! Number of active cycles completed so far
    unsigned int d_cycle;

    //! Store individual keff estimators over all (active + inactive) cycles
    Vec_Dbl d_all_keff;

    //! Accumulated first moment of keff for calculating average
    double d_keff_sum;

    //! Accumulated second moment of keff for calculating variance
    double d_keff_sum_sq;

  public:
    // Kcode solver should construct this with initial keff estimate
    Keff_Tally(double keff_init, SP_Physics physics);

    // >>> ACCESSORS

    //! Access all keff estimators, both active and inactive
    const Vec_Dbl& all_keff() const { return d_all_keff; }

    //! Obtain keff estimate from this cycle
    double latest() const { return d_keff_cycle; }

    // Calculate average keff over active cycles
    double mean() const;

    // Calculate variance of average keff estimate over active cycles
    double variance() const;

    //! Number of cycles since resetting
    unsigned int cycle_count() const { return d_cycle; }

    // >>> ACCESSORS FOR TESTING

    //! Obtain first moment of keff (for testing purposes)
    double keff_sum() const { return d_keff_sum; }

    //! Obtain second moment of keff (for testing purposes)
    double keff_sum_sq() const { return d_keff_sum_sq; }

    // >>> DERIVED INTERFACE

    // Track particle and do tallying.
    virtual void accumulate(double step, const Particle_t &p) override final;

    //! Accumulate first and second moments
    virtual void end_history() final { /* * */ }

    //! Do post-processing on first and second moments
    virtual void finalize(double num_histories) final { /* * */ }

    // Begin active cycles in a kcode calculation.
    virtual void begin_active_cycles() override final;

    // Begin a new cycle in a kcode calculation.
    virtual void begin_cycle() override final;

    // End a cycle in a kcode calculation.
    virtual void end_cycle(double num_particles) override final;

    // Clear/re-initialize all tally values between solves.
    virtual void reset() final;
};

} // end namespace profugus

#endif // mc_Keff_Tally_hh

//---------------------------------------------------------------------------//
//                 end of Keff_Tally.hh
//---------------------------------------------------------------------------//
