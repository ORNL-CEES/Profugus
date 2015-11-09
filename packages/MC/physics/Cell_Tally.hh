//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Cell_Tally.hh
 * \author Thomas M. Evans
 * \date   Fri Aug 21 15:29:14 2015
 * \brief  Cell_Tally class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Cell_Tally_hh
#define MC_mc_Cell_Tally_hh

#include <unordered_map>
#include <vector>
#include <memory>
#include <utility>

#include "Cell_Tally_State.hh"
#include "Tally.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Cell_Tally
 * \brief Do pathlength cell tallies.
 */
/*!
 * \example mc/test/tstCell_Tally.cc
 *
 * Test of Cell_Tally.
 */
//===========================================================================//

class Cell_Tally : public Pathlength_Tally
{
  private:

    // Base typedef.
    typedef Pathlength_Tally Base;

    // Atomic tally moment.
    class Atomic_Moment
    {
      private:
	double d_first_moment;
	double d_second_moment;
      public:
	Atomic_Moment() { /* ... */ }
	Atomic_Moment( const double first_moment, const double second_moment )
	    : d_first_moment( first_moment )
	    , d_second_moment( second_moment )
	{ /* ... */ }
	Atomic_Moment( const Atomic_Moment& rhs )
	    : d_first_moment( rhs.d_first_moment )
	    , d_second_moment( rhs.d_second_moment )
	{ /* ... */ }
	inline double first() const { return d_first_moment; }
	inline double& first() { return d_first_moment; }
	inline double second() const { return d_second_moment; }
	inline double& second() { return d_second_moment; }
	inline void atomic_sum_into( const double first, const double second )
	{
	    d_first_moment += first;
	    d_second_moment += second;
	}
    };

  public:
    //@{
    //! Typedefs.
    typedef Cell_Tally_State                 State_t;
    typedef Physics_t::Geometry_t            Geometry_t;
    typedef std::shared_ptr<Geometry_t>      SP_Geometry;
    typedef Atomic_Moment                    Moments;
    typedef std::unordered_map<int, Moments> Result;
    //@}

  private:
    // >>> DATA

    // Geometry.
    SP_Geometry d_geometry;

    // Map of tally moments.
    Result d_tally;

  public:
    // Constructor.
    Cell_Tally(SP_Physics physics);

    // Add tally cells.
    void set_cells(const std::vector<int> &cells);

    // Get tally results.
    const Result& results() const { return d_tally; }

    // >>> TALLY INTERFACE

    // Accumulate first and second moments
    void end_history(const Particle_t &p);

    // Do post-processing on first and second moments
    void finalize(double num_particles);

    // Clear/re-initialize all tally values between solves
    void reset();

    // Output results.
    void output(const std::string &out_file);

    // >>> PATHLENGTH TALLY INTERFACE

    // Track particle and tally.
    void accumulate( const double step, Particle_t &p) const;

  private:
    // >>> IMPLEMENTATION

    // Index in the particles tally state.
    const unsigned int d_state_idx;

    // Clear local values.
    void clear_local();
};

//---------------------------------------------------------------------------//

} // end namespace profugus

#endif // MC_mc_Cell_Tally_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Cell_Tally.hh
//---------------------------------------------------------------------------//
