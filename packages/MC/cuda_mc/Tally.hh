//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Tally.hh
 * \author Thomas M. Evans
 * \date   Wed May 14 15:10:09 2014
 * \brief  Tally base class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Tally_hh
#define cuda_mc_Tally_hh

#include <memory>
#include <string>

#include "cuda_utils/Shared_Device_Ptr.hh"

#include "Particle_Vector.hh"
#include "Definitions.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Tally
 * \brief Base class for tallies.
 */
//===========================================================================//

template <class Geometry>
class Tally
{
  public:
    // Destructor.
    virtual ~Tally() { /* ... */ }

    // >>> PUBLIC INTERFACE

    //! Query if this tally is on during inactive cycles.
    virtual bool inactive_cycle_tally() const { return true; }

    //! Accumulate first and second moments
    virtual void end_history() { /* * */ }

    //! Do post-processing on first and second moments
    virtual void finalize(double num_particles) { /* * */ }

    //! Begin active cycles in a kcode calculation (default no-op)
    virtual void begin_active_cycles() { /* * */ }

    //! Begin a new cycle in a kcode calculation (default no-op)
    virtual void begin_cycle() { /* * */ }

    //! End a cycle in a kcode calculation (default no-op)
    virtual void end_cycle(double num_particles) { /* * */ }

    //! Clear/re-initialize all tally values between solves
    virtual void reset() { /* * */ }
};

//---------------------------------------------------------------------------//
/*!
 * \class Source_Tally
 * \brief Defines source tally interfaces.
 */
template <class Geometry>
class Source_Tally : public Tally<Geometry>
{
    typedef Tally<Geometry>                 Base;

  public:

    // Destructor.
    virtual ~Source_Tally() { /* ... */ }

    // >>> TALLY INTERFACE

    //! Tally events at particle birth.
    virtual void birth(
	const cuda::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles ) = 0;
};

//---------------------------------------------------------------------------//
/*!
 * \class Pathlength_Tally
 * \brief Defines source tally interfaces.
 */
template <class Geometry>
class Pathlength_Tally : public Tally<Geometry>
{
    typedef Tally<Geometry>                 Base;

  public:

    // Destructor.
    virtual ~Pathlength_Tally() { /* ... */ }

    // >>> TALLY INTERFACE

    //! Track particle and tally.
    virtual void accumulate(
	const cuda::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles 
	) = 0;
};

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Tally_hh

//---------------------------------------------------------------------------//
//                 end of Tally.hh
//---------------------------------------------------------------------------//
