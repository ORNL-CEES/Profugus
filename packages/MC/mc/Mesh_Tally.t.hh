//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Mesh_Tally.t.hh
 * \author Thomas M. Evans
 * \date   Fri Aug 21 15:29:14 2015
 * \brief  Mesh_Tally class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Mesh_Tally_t_hh
#define MC_mc_Mesh_Tally_t_hh

#include <cmath>
#include <algorithm>

#include "Utils/comm/global.hh"
#include "Utils/harness/Soft_Equivalence.hh"
#include "Mesh_Tally.hh"
#include "utils/Serial_HDF5_Writer.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Mesh_Tally<Geometry>::Mesh_Tally(RCP_Std_DB db, SP_Physics physics)
    : Base(physics, false)
    , d_geometry(b_physics->get_geometry())
    , d_db(db)
{
    REQUIRE(d_geometry);

    // set the tally name
    set_name("mesh");

    REQUIRE( d_db->isType<std::string>("problem_name") );
    d_filename = d_db->get<std::string>("problem_name") + "_mesh_tally.h5";

    // reset tally
    reset();

    d_cycle_output = d_db->get("do_cycle_output",false);
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Add a list of cells to tally.
 */
template <class Geometry>
void Mesh_Tally<Geometry>::set_mesh(SP_Mesh mesh)
{
    REQUIRE( mesh );

    d_mesh = mesh;

    // Resize result vector
    d_tally.resize(mesh->num_cells(),{0.0,0.0});
    d_hist.resize(mesh->num_cells(),0.0);
    d_cycle_tally.resize(mesh->num_cells(),0.0);

#ifdef USE_HDF5
    Serial_HDF5_Writer writer;
    writer.open(d_filename);
    writer.write("x_edges",d_mesh->edges(def::I));
    writer.write("y_edges",d_mesh->edges(def::J));
    writer.write("z_edges",d_mesh->edges(def::K));
    writer.close();
#endif
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*
 * \brief Accumulate first and second moments.
 */
template <class Geometry>
void Mesh_Tally<Geometry>::end_history()
{
    REQUIRE( d_tally.size() == d_mesh->num_cells() );
    REQUIRE( d_hist.size()  == d_mesh->num_cells() );

    // Iterate through local tally results and add them to the permanent
    // results
    for (int cell = 0; cell < d_mesh->num_cells(); ++cell )
    {
        d_tally[cell].first  += d_hist[cell];
        d_tally[cell].second += d_hist[cell]*d_hist[cell];
        d_cycle_tally[cell] += d_hist[cell];
    }

    // Clear the local tally
    clear_local();
}

//---------------------------------------------------------------------------//
/*
 * \brief Begin new cycle
 */
template <class Geometry>
void Mesh_Tally<Geometry>::begin_cycle()
{
    REQUIRE( d_cycle_tally.size() == d_mesh->num_cells() );

    std::fill(d_cycle_tally.begin(),d_cycle_tally.end(),0.0);
}

//---------------------------------------------------------------------------//
/*
 * \brief End cycle
 */
template <class Geometry>
void Mesh_Tally<Geometry>::end_cycle(double num_particles)
{
    if (d_cycle_output)
    {
        REQUIRE( d_cycle_tally.size() == d_mesh->num_cells() );
        REQUIRE( num_particles > 0.0 );

        profugus::global_sum(&d_cycle_tally[0], d_cycle_tally.size());

        for (int cell = 0; cell < d_cycle_tally.size(); ++cell)
        {
            CHECK(cell<d_mesh->num_cells());
            double nrm_factor = num_particles * d_mesh->volume(cell);
            CHECK(nrm_factor>0.0);
            d_cycle_tally[cell] /= nrm_factor;
        }

#ifdef USE_HDF5
        Serial_HDF5_Writer writer;
        writer.open(d_filename,HDF5_IO::APPEND);

        std::ostringstream m;
        m << "cycle_" << d_cycle;
        writer.begin_group(m.str());

        writer.write("cycle_flux",d_cycle_tally);

        writer.end_group();
        writer.close();
#endif
    }

    d_cycle++;
}

//---------------------------------------------------------------------------//
/*
 * \brief Do post-processing on first and second moments.
 */
template <class Geometry>
void Mesh_Tally<Geometry>::finalize(double num_particles)
{
    REQUIRE(num_particles > 1);

    // Do a global reduction on moments
    std::vector<double> first( d_tally.size(),  0.0);
    std::vector<double> second(d_tally.size(), 0.0);

    // Write tally results into separate arrays for reduction
    for( int cell = 0; cell < d_tally.size(); ++cell )
    {
        first[cell]  = d_tally[cell].first;
        second[cell] = d_tally[cell].second;
    }

    // Do global reductions on the moments
    profugus::global_sum(&first[0],  d_tally.size());
    profugus::global_sum(&second[0], d_tally.size());

    // Store 1/N
    double inv_N = 1.0 / num_particles;

    for( int cell = 0; cell < d_tally.size(); ++cell )
    {
        double inv_V = 1.0 / d_mesh->volume(cell);
        CHECK( inv_V > 0.0 );

        // Calculate means for this cell
        double avg_l  = first[cell]  * inv_N;
        double avg_l2 = second[cell] * inv_N;

        // Get a reference to the moments
        auto &moments = d_tally[cell];

        // Store the sample mean
        moments.first = avg_l * inv_V;

        // Calculate the variance
        double var = (avg_l2 - avg_l * avg_l) / (num_particles - 1);

        // Store the error of the sample mean
        moments.second = std::sqrt(var) * inv_V;

        // Store individual vectors for output purposes
        first[cell]  = moments.first;
        second[cell] = moments.second;
    }

#ifdef USE_HDF5
    Serial_HDF5_Writer writer;
    writer.open(d_filename,HDF5_IO::APPEND);
    writer.write("flux_mean",first);
    writer.write("flux_std_dev",second);
    writer.close();
#endif
}

//---------------------------------------------------------------------------//
/*
 * \brief Clear/re-initialize all tally values between solves.
 */
template <class Geometry>
void Mesh_Tally<Geometry>::reset()
{
    // Clear the local tally
    clear_local();

    d_cycle = 0;
}

//---------------------------------------------------------------------------//
/*
 * \brief Track particle and tally..
 */
template <class Geometry>
void Mesh_Tally<Geometry>::accumulate(double            step,
                                         const Particle_t &p)
{
    using def::I; using def::J; using def::K;

    // Get the cell index
    Mesh::Dim_Vector ijk;
    const Mesh::Space_Vector &xyz = p.geo_state().d_r;
    d_mesh->find_upper(xyz,ijk);

    // Check boundaries for particles pushed just outside domain
    using profugus::soft_equiv;
    constexpr double tol = 1e-12;
    for (auto dir : {I, J, K})
    {
        if (soft_equiv(xyz[dir],d_mesh->low_corner(dir),tol))
            ijk[dir] = 0;
        else if (soft_equiv(xyz[dir],d_mesh->high_corner(dir),tol))
            ijk[dir] = d_mesh->num_cells_along(dir)-1;
    }

    Mesh::size_type cell;
    bool found = d_mesh->index( ijk[I], ijk[J], ijk[K], cell );
    if( found )
    {
        REQUIRE( cell >= 0 && cell < d_mesh->num_cells() );

        // Tally for the history
        d_hist[cell] += p.wt() * step;
    }
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*
 * \brief Clear local values.
 */
template <class Geometry>
void Mesh_Tally<Geometry>::clear_local()
{
    // Clear the local tally
    std::fill(d_hist.begin(), d_hist.end(), 0.0);
}

} // end namespace profugus

#endif // MC_mc_Mesh_Tally_t_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Mesh_Tally.t.hh
//---------------------------------------------------------------------------//
