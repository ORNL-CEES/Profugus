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

#include <vector>
#include <memory>

#include "Tally.hh"
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
    typedef Physics                             Physics_t;
    typedef Physics_t::Geometry_t               Geometry_t;
    typedef Physics_t::Particle_t               Particle_t;
    typedef Tally                               Tally_t;
    typedef Pathlength_Tally                    Pathlength_Tally_t;
    typedef Source_Tally                        Source_Tally_t;
    typedef Compound_Tally                      Compound_Tally_t;
    typedef std::shared_ptr<Tally_t>            SP_Tally;
    typedef std::shared_ptr<Pathlength_Tally_t> SP_Pathlength_Tally;
    typedef std::shared_ptr<Source_Tally_t>     SP_Source_Tally;
    typedef std::shared_ptr<Compound_Tally_t>   SP_Compound_Tally;
    typedef std::shared_ptr<Geometry_t>         SP_Geometry;
    typedef std::shared_ptr<Physics_t>          SP_Physics;
    typedef std::vector<SP_Tally>               Vec_Tallies;
    //@}

  private:
    // >>> DATA

    // Geometry and physics.
    SP_Geometry d_geometry;
    SP_Physics  d_physics;

    // Vector of all tallies (assembled during build).
    Vec_Tallies d_tallies;

    // Persistent source and pathlength tallies.
    std::vector<SP_Pathlength_Tally> d_pl;
    std::vector<SP_Source_Tally>     d_src;
    std::vector<SP_Compound_Tally>   d_comp;

  public:
    // Constructor.
    Tallier();

    // Set the geometry and physics classes.
    void set(SP_Geometry geometry, SP_Physics physics);

    //! Get the geometry we use.
    SP_Geometry geometry() const { return d_geometry; }

    //! Get the physics we use.
    SP_Physics physics() const { return d_physics; }

    // Add tallies.
    void add_pathlength_tally(SP_Pathlength_Tally tally);
    void add_source_tally(SP_Source_Tally tally);
    void add_compound_tally(SP_Compound_Tally tally);

    //@{
    //! Number of tallies.
    auto num_tallies() const -> decltype(d_tallies.size())
    {
        return d_tallies.size();
    }
    auto num_pathlength_tallies() const -> decltype(d_pl.size())
    {
        return d_pl.size();
    }
    auto num_source_tallies() const -> decltype(d_src.size())
    {
        return d_src.size();
    }
    auto num_compound_tallies() const -> decltype(d_comp.size())
    {
        return d_comp.size();
    }
    //@}

    //@{
    //! Iterate through all tallies.
    auto begin() -> decltype(d_tallies.begin()) { return d_tallies.begin(); }
    auto end()   -> decltype(d_tallies.end())   { return d_tallies.end(); }
    //@}

    // Initialize internal data structures after adding tallies.
    void build();

    // >>> TALLY OPERATIONS

    // Process path-length tally events.
    void path_length(double step, const Particle_t &p);

    // Tally any source events.
    void source(const Particle_t &p);

    // Tell the tallies to begin active kcode cycles
    void begin_active_cycles();

    // Tell the tallies to begin a new cycle in a kcode calculation
    void begin_cycle();

    // Tell the tallies to end a cycle in a kcode calculation
    void end_cycle(double num_particles);

    // Perform all end-history tally tasks.
    void end_history();

    // Finalize tallies.
    void finalize(double num_particles);

    // Reset tallies.
    void reset();

    // Swap two talliers.
    void swap(Tallier &rhs);

    // >>> ACCESSORS

    //! Whether we've called "build" with the current number of tallies
    bool is_built() const { return d_build_phase == BUILT; }

    //! Whether we've called "finalize" with the given tallies
    bool is_finalized() const { return d_build_phase == FINALIZED; }

  private:
    // IMPLEMENTATION

    // Prune tallies for doubles.
    template<class Vec_T>
    void prune(Vec_T &tallies);

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

//---------------------------------------------------------------------------//
// NON-MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Swap two talliers.
 *
 * This is useful for deactivating tallying (like in a KCODE problem).
 *
 * This provides a std-like swap solution using Koenig namespace lookup.
 */
inline void swap(Tallier &a,
                 Tallier &b)
{
    a.swap(b);
}

} // end namespace profugus

#endif // mc_Tallier_hh

//---------------------------------------------------------------------------//
//                 end of Tallier.hh
//---------------------------------------------------------------------------//
