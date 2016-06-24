//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Mesh_Tally.hh
 * \author Thomas M. Evans
 * \date   Fri Aug 21 15:29:14 2015
 * \brief  Mesh_Tally class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Mesh_Tally_hh
#define MC_mc_Mesh_Tally_hh

#include <vector>
#include <memory>
#include <utility>

#include "Tally.hh"
#include "geometry/Cartesian_Mesh.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Mesh_Tally
 * \brief Do pathlength mesh tallies on fission source.
 */
/*!
 * \example mc/test/tstMesh_Tally.cc
 *
 * Test of Mesh_Tally.
 */
//===========================================================================//

template <class Geometry>
class Mesh_Tally : public Pathlength_Tally<Geometry>
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
    typedef std::vector<Moments>             Result;
    typedef std::shared_ptr<Physics_t>       SP_Physics;
    typedef profugus::Cartesian_Mesh         Mesh;
    typedef std::shared_ptr<Mesh>            SP_Mesh;
    typedef Teuchos::ParameterList           ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t>    RCP_Std_DB;
    //@}

  private:
    // >>> DATA

    using Base::b_physics;

    // Geometry.
    SP_Geometry d_geometry;

    // Cartesian mesh
    SP_Mesh d_mesh;

    // Map of tally moments.
    Result d_tally;

    // Should fluxes be written every cycle?
    bool d_cycle_output;

    // Cycle tally
    std::vector<double> d_cycle_tally;
    int d_cycle;

    // Parameters
    RCP_Std_DB d_db;

    // HDF5 output filename
    std::string d_filename;

  public:
    // Constructor.
    Mesh_Tally(RCP_Std_DB db, SP_Physics physics);

    // Add tally mesh.
    void set_mesh(SP_Mesh mesh);

    // Get tally results.
    const Result& results() const { return d_tally; }

    // >>> TALLY INTERFACE

    // Accumulate first and second moments
    void end_history();

    // Begin new cycle
    void begin_cycle();

    // End cycle
    void end_cycle(double num_particles);

    // Do post-processing on first and second moments
    void finalize(double num_particles);

    // Clear/re-initialize all tally values between solves
    void reset();

    // >>> PATHLENGTH TALLY INTERFACE

    // Track particle and tally.
    void accumulate(double step, const Particle_t &p);

  private:
    // >>> IMPLEMENTATION

    typedef std::vector<double> History_Tally;

    // Clear local values.
    void clear_local();

    // Tally for a history.
    History_Tally d_hist;
};

//---------------------------------------------------------------------------//

} // end namespace profugus

#endif // MC_mc_Mesh_Tally_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Mesh_Tally.hh
//---------------------------------------------------------------------------//
