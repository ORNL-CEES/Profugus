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

template <class Geometry>
class Cell_Tally : public Pathlength_Tally<Geometry>
{
    typedef Pathlength_Tally<Geometry> Base;
    using Base::set_name;

  public:
    //@{
    //! Typedefs.
    typedef Geometry                         Geometry_t;
    typedef Physics<Geometry>                Physics_t;
    typedef typename Physics_t::Particle_t   Particle_t;
    typedef std::shared_ptr<Geometry_t>      SP_Geometry;
    typedef std::pair<double, double>        Moments;
    typedef std::unordered_map<int, Moments> Result;
    typedef std::shared_ptr<Physics_t>       SP_Physics;
    typedef Teuchos::ParameterList           ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t>    RCP_Std_DB;
    //@}

  private:
    // >>> DATA

    using Base::b_physics;

    // Geometry.
    SP_Geometry d_geometry;

    // Map of tally moments.
    Result d_tally;

    // Database
    RCP_Std_DB d_db;

  public:
    // Constructor.
    Cell_Tally(RCP_Std_DB db, SP_Physics physics);

    // Add tally cells.
    void set_cells(const std::vector<int> &cells);

    // Get tally results.
    const Result& results() const { return d_tally; }

    // >>> TALLY INTERFACE

    // Begin new cycle
    void begin_cycle();

    // End cycle
    void end_cycle(double num_particles);

    // Accumulate first and second moments
    void end_history();

    // Do post-processing on first and second moments
    void finalize(double num_particles);

    // Clear/re-initialize all tally values between solves
    void reset();

    // >>> PATHLENGTH TALLY INTERFACE

    // Track particle and tally.
    void accumulate(double step, const Particle_t &p);

  private:
    // >>> IMPLEMENTATION

    typedef std::unordered_map<int, double> History_Tally;

    // Clear local values.
    void clear_local();

    // Cycle counter
    int d_cycle;

    // Map for reordering fluxes from unordered map
    std::vector<int> d_cell_map;

    // Filename for HDF5 output
    std::string d_outfile;

    // Tally for a history.
    History_Tally d_hist;

    // Should we write fluxes at every cycle
    bool d_cycle_output;

    // Tally for single cycle
    History_Tally d_cycle_tally;
};

//---------------------------------------------------------------------------//

} // end namespace profugus

#endif // MC_mc_Cell_Tally_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Cell_Tally.hh
//---------------------------------------------------------------------------//
