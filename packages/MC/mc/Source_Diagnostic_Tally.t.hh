//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/Source_Diagnostic_Tally.t.hh
 * \author Thomas M. Evans
 * \date   Tue Dec 09 16:32:04 2014
 * \brief  Source_Diagnostic_Tally member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Source_Diagnostic_Tally_t_hh
#define MC_mc_Source_Diagnostic_Tally_t_hh

#include <algorithm>
#include <sstream>
#include <cmath>

#include "harness/DBC.hh"
#include "harness/Warnings.hh"
#include "comm/global.hh"
#include "utils/Definitions.hh"
#include "Source_Diagnostic_Tally.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Source_Diagnostic_Tally<Geometry>::Source_Diagnostic_Tally(
        RCP_Std_DB       db,
        SP_Physics       physics,
        SP_Mesh_Geometry mesh,
        bool             inactive)
    : Base(physics, inactive)
    , d_mesh(mesh)
    , d_geometry(b_physics->get_geometry())
    , d_source_density(d_mesh->num_cells(), 0.0)
    , d_cycle_ctr(0)
    , d_num_per_cycle(0)
{
    REQUIRE(d_mesh);
    REQUIRE(!db.is_null());
    REQUIRE(d_geometry);
    REQUIRE(db->isParameter("problem_name"));
    REQUIRE(db->isParameter("num_cycles"));

    // set the name
    Base::set_name("fission_source");

#ifndef USE_HDF5
    ADD_WARNING("HDF5 not available in this build, turning source diagnostic"
                << " tally off");
#else

    // make the source diagnostic output file
    d_filename = db->get<std::string>("problem_name") + "_fs.h5";

    // make the initial file
    d_writer.open(d_filename);

    // write the mesh
    d_writer.begin_group("mesh");

    // get the underlying cartesian mesh
    const auto &cart_mesh = d_mesh->mesh();

    // write the edges
    d_writer.write("x", cart_mesh.edges(def::X));
    d_writer.write("y", cart_mesh.edges(def::Y));
    d_writer.write("z", cart_mesh.edges(def::Z));

    d_writer.end_group();

    d_writer.close();

#endif
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Tally particles at birth.
 */
template <class Geometry>
void Source_Diagnostic_Tally<Geometry>::birth(const Particle_t &p)
{
#ifdef USE_HDF5
    // get the particle's geometric state
    const auto &geo_state = p.geo_state();

    // initialize the particle in the fission-matrix mesh
    d_mesh->initialize(d_geometry->position(geo_state),
                       d_geometry->direction(geo_state),
                       d_state);

    // determine the birth cell in the fission matrix mesh
    auto mesh_idx = d_mesh->cell(d_state);
    CHECK(mesh_idx >= 0 && mesh_idx < d_mesh->num_cells());

    // tally the source (particle) density
    d_source_density[mesh_idx] += p.wt();

    // count up the particles in this cycle
    ++d_num_per_cycle;
#endif
}

//---------------------------------------------------------------------------//
/*!
 * \brief End a cycle in a kcode calculation.
 */
template <class Geometry>
void Source_Diagnostic_Tally<Geometry>::end_cycle(double num_particles)
{
#ifdef USE_HDF5
    REQUIRE(d_source_density.size() == d_mesh->num_cells());

    // reduce the tally across all sets
    profugus::global_sum(d_source_density.data(), d_source_density.size());

    // get the underlying cartesian mesh
    const auto &cart_mesh = d_mesh->mesh();
    CHECK(cart_mesh.num_cells() == d_mesh->num_cells());

    // normalize so that s\dot s = 1
    double norm = 0.0;
    for (int cell = 0; cell < d_source_density.size(); ++cell)
    {
        d_source_density[cell] /= cart_mesh.volume(cell);
        norm += d_source_density[cell] * d_source_density[cell];
    }

    // calculate the normalization
    norm = 1.0 / std::sqrt(norm);

    // apply the normalization
    for (auto &s : d_source_density)
    {
        s *= norm;
    }

    // open the file - writing is only on proc 0
    d_writer.open(d_filename, HDF5_IO::APPEND, 0);

    // make the cycle group
    std::ostringstream m;
    m << "cycle_" << d_cycle_ctr;
    d_writer.begin_group(m.str());

    // write the number of particles
    profugus::global_sum(&d_num_per_cycle, 1);
    d_writer.write("num_particles", d_num_per_cycle);

    // write the normalized source
    d_writer.write("source_density", d_source_density);

    // end the group and close
    d_writer.end_group();
    d_writer.close();

    // reset the source tally for the next cycle
    reset();

    // update the cycle counter
    ++d_cycle_ctr;

#endif
}

//---------------------------------------------------------------------------//
/*!
 * \brief Clear/re-initialize all tally values.
 */
template <class Geometry>
void Source_Diagnostic_Tally<Geometry>::reset()
{
    std::fill(d_source_density.begin(), d_source_density.end(), 0.0);
    d_num_per_cycle = 0;
}

} // end namespace profugus

#endif // MC_mc_Source_Diagnostic_Tally_t_hh

//---------------------------------------------------------------------------//
//                 end of Source_Diagnostic_Tally.t.hh
//---------------------------------------------------------------------------//
