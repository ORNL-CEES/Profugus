//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Tallier.hh
 * \author Thomas M. Evans and Seth Johnson
 * \date   Mon May 12 12:15:30 2014
 * \brief  Tallier class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Tallier_hh
#define mc_Tallier_hh

#include <memory>

#include "Physics.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Tallier
 * \brief Do tally operations.
 */
/*!
 * \example mc/test/tstTallier.cc
 *
 * Test of Tallier.
 */
//===========================================================================//

class Tallier
{
  public:
    //@{
    //! Typedefs.
    typedef Physics                     Physics_t;
    typedef Physics_t::Geometry_t       Geometry_t;
    typedef Physics_t::Particle_t       Particle_t;
    typedef std::shared_ptr<Geometry_t> SP_Geometry;
    typedef std::shared_ptr<Physics_t>  SP_Physics;
    //@}

  private:
    // >>> DATA

    // Geometry and physics.
    SP_Geometry d_geometry;
    SP_Physics  d_physics;

  public:
    // Constructor.
    Tallier();

    // Set the geometry and physics classes.
    void set(SP_Geometry geometry, SP_Physics physics);

    // Initialize internal data structures after adding tallies.
    void build();

    // Process path-length tally events.
    void path_length(double step, const Particle_t &p);

    // Tally any source events.
    void source(const Particle_t &p);

    // Perform all end-history tally tasks.
    void end_history();

    // Finalize tallies.
    void finalize(int Np);

    // Reset tallies.
    void reset();

    // >>> ACCESSORS

    //! Whether we've called "build" with the current number of tallies
    bool is_built() const { return d_build_phase == BUILT; }

    //! Whether we've called "finalize" with the given tallies
    bool is_finalized() const { return d_build_phase == FINALIZED; }

  private:
    // IMPLEMENTATION

    //! Phases of construction, for error checking
    enum Build_Phase
    {
        CONSTRUCTED = 0,//!< after construction is complete
        ASSIGNED,       //!< after assigning geometry and physics
        BUILT,          //!< after the call to build()
        FINALIZED       //!< after the call to finalize()
    };

    // Build phase.
    Build_Phase d_build_phase;
};

} // end namespace profugus

#endif // mc_Tallier_hh

//---------------------------------------------------------------------------//
//                 end of Tallier.hh
//---------------------------------------------------------------------------//
