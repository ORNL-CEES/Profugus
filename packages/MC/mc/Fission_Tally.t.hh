//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Fission_Tally.t.hh
 * \author Thomas M. Evans
 * \date   Fri Aug 21 15:29:14 2015
 * \brief  Fission_Tally class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Fission_Tally_t_hh
#define MC_mc_Fission_Tally_t_hh

#include <cmath>
#include <algorithm>

#include "Utils/comm/global.hh"
#include "Fission_Tally.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Fission_Tally<Geometry>::Fission_Tally(SP_Physics physics)
    : Base(physics, false)
    , d_geometry(b_physics->get_geometry())
{
    REQUIRE(d_geometry);

    // set the tally name
    set_name("fission");

    // reset tally
    reset();
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Add a list of cells to tally.
 */
template <class Geometry>
void Fission_Tally<Geometry>::set_mesh(SP_Mesh mesh)
{
    REQUIRE( mesh );

    d_mesh = mesh;

    // Resize result vector
    d_tally.resize(mesh->num_cells(),{0.0,0.0});
    d_hist.resize( mesh->num_cells(),0.0);
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*
 * \brief Accumulate first and second moments.
 */
template <class Geometry>
void Fission_Tally<Geometry>::end_history()
{
    REQUIRE( d_tally.size() == d_mesh->num_cells() );
    REQUIRE( d_hist.size()  == d_mesh->num_cells() );

    // Iterate through local tally results and add them to the permanent
    // results
    for (int cell = 0; cell < d_mesh->num_cells(); ++cell )
    {
        d_tally[cell].first  += d_hist[cell];
        d_tally[cell].second += d_hist[cell]*d_hist[cell];
    }

    // Clear the local tally
    clear_local();
}

//---------------------------------------------------------------------------//
/*
 * \brief Do post-processing on first and second moments.
 */
template <class Geometry>
void Fission_Tally<Geometry>::finalize(double num_particles)
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
    double inv_N = 1.0 / static_cast<double>(num_particles);

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
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Clear/re-initialize all tally values between solves.
 */
template <class Geometry>
void Fission_Tally<Geometry>::reset()
{
    // Clear the local tally
    clear_local();
}

//---------------------------------------------------------------------------//
/*
 * \brief Track particle and tally..
 */
template <class Geometry>
void Fission_Tally<Geometry>::accumulate(double            step,
                                         const Particle_t &p)
{
    // Get the cell index
    Mesh::Dim_Vector ijk;
    const Mesh::Space_Vector &xyz = p.geo_state().d_r;
    d_mesh->find(xyz,ijk);
    int cell = d_mesh->index( ijk[def::I], ijk[def::J], ijk[def::K] );

    // Tally for the history
    d_hist[cell] += p.wt() * step * b_physics->total(physics::NU_FISSION,p);
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*
 * \brief Clear local values.
 */
template <class Geometry>
void Fission_Tally<Geometry>::clear_local()
{
    // Clear the local tally
    std::fill(d_hist.begin(), d_hist.end(), 0.0);
}

} // end namespace profugus

#endif // MC_mc_Fission_Tally_t_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Fission_Tally.t.hh
//---------------------------------------------------------------------------//
