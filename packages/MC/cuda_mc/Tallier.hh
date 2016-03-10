//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Tallier.hh
 * \author Stuart Slattery
 * \date   Mon May 12 12:15:30 2014
 * \brief  Tallier class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Tallier_hh
#define cuda_mc_Tallier_hh

#include <vector>
#include <memory>

#include "cuda_utils/Shared_Device_Ptr.hh"

#include "Tally.hh"
#include "Particle_Vector.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Tallier
 * \brief Do tally operations.
 */
/*!
 * \example cuda_mc/test/tstTallier.cc
 *
 * Test of Tallier.
 */
//===========================================================================//

template <class Geometry>
class Tallier
{
  public:
    //@{
    //! Typedefs.
    typedef Geometry                            Geometry_t;
    typedef Tally<Geometry_t>                   Tally_t;
    typedef Pathlength_Tally<Geometry_t>        Pathlength_Tally_t;
    typedef Source_Tally<Geometry_t>            Source_Tally_t;
    typedef std::shared_ptr<Tally_t>            SP_Tally;
    typedef std::shared_ptr<Pathlength_Tally_t> SP_Pathlength_Tally;
    typedef std::shared_ptr<Source_Tally_t>     SP_Source_Tally;
    typedef std::vector<SP_Tally>               Vec_Tallies;
    typedef Particle_Vector<Geometry>           Particle_Vector_t;
    //@}

  private:
    // >>> DATA

    // Vector of all tallies (assembled during build).
    Vec_Tallies d_tallies;

    // Persistent source and pathlength tallies.
    std::vector<SP_Pathlength_Tally> d_pl;
    std::vector<SP_Source_Tally>     d_src;

  public:
    // Constructor.
    Tallier();

    // Add tallies.
    void add_pathlength_tally(SP_Pathlength_Tally tally);
    void add_source_tally(SP_Source_Tally tally);

    //@{
    //! Number of tallies.
    int num_tallies() const { return d_tallies.size(); }
    int num_pathlength_tallies() const { return d_pl.size(); }
    int num_source_tallies() const { return d_src.size(); }
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
    void path_length(
	const cuda::Shared_Device_Ptr<Particle_Vector_t>& particles );

    // Tally any source events.
    void source(
	const cuda::Shared_Device_Ptr<Particle_Vector_t>& particles );

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

    //! Phases of construction, for error checking
    enum Build_Phase
    {
        CONSTRUCTED = 0,//!< after construction is complete
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
template <class Geometry>
inline void swap(Tallier<Geometry> &a,
                 Tallier<Geometry> &b)
{
    a.swap(b);
}

} // end namespace cuda_profugus

#endif // cuda_mc_Tallier_hh

//---------------------------------------------------------------------------//
//                 end of Tallier.hh
//---------------------------------------------------------------------------//
