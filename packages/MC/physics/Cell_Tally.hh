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
    typedef Pathlength_Tally Base;

  public:
    //@{
    //! Typedefs.
    typedef Cell_Tally_State                 State_t;
    typedef Physics_t::Geometry_t            Geometry_t;
    typedef std::shared_ptr<Geometry_t>      SP_Geometry;
    typedef std::pair<double, double>        Moments;
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
    void accumulate(double step, const Particle_t &p);

  private:
    // >>> IMPLEMENTATION

    // Index in the particles tally state.
    const unsigned int d_state_idx;

    typedef std::unordered_map<int, double> History_Tally;

    // Clear local values.
    void clear_local();

    // Tally for a history.
    History_Tally d_hist;
};

//---------------------------------------------------------------------------//

} // end namespace profugus

#endif // MC_mc_Cell_Tally_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Cell_Tally.hh
//---------------------------------------------------------------------------//
