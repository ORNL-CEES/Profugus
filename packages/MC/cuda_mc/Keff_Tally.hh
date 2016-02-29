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

#include "Physics.hh"
#include "Particle_Vector.hh"

#include "utils/Definitions.hh"

#include "cuda_utils/Shared_Device_Ptr.hh"

namespace cuda_profugus
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

template <class Geometry>
class Keff_Tally
{
  public:
    //@{
    //! Misc typedefs.
    typedef Particle_Vector<Geometry> Particle_Vector_t;
    typedef Physics<Geometry> Physics_t;
    typedef def::Vec_Dbl Vec_Dbl;
    //@}

  private:
    // >>> DATA

    //! Physics.
    cuda::Shared_Device_Ptr<Physics_t> d_physics;
    
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

    //! keff host work vector.
    std::vector<double> d_keff_host;

    //! keff work vector on-device.
    double* d_keff_device;

  public:
    
    // Kcode solver should construct this with initial keff estimate
    Keff_Tally( const double keff_init,
		const cuda::Shared_Device_Ptr<Physics_t>& physics );

    // Destructor.
    ~Keff_Tally();
    
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
    void accumulate( cuda::Shared_Device_Ptr<Particle_Vector_t>& particles );

    // Begin active cycles in a kcode calculation.
    void begin_active_cycles();

    // Begin a new cycle in a kcode calculation.
    void begin_cycle();

    // End a cycle in a kcode calculation.
    void end_cycle(double num_particles);

    // Clear/re-initialize all tally values between solves.
    void reset();

    // >>> SETTERS

    //! Set the latest keff in the tally.
    void set_keff(double k) { d_keff_cycle = k; }
};

} // end namespace cuda_profugus

#endif // mc_Keff_Tally_hh

//---------------------------------------------------------------------------//
//                 end of Keff_Tally.hh
//---------------------------------------------------------------------------//
