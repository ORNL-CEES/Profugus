//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Cell_Tally.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 21 15:29:14 2015
 * \brief  Cell_Tally class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <algorithm>

#include "comm/global.hh"
#include "utils/Serial_HDF5_Writer.hh"
#include "Cell_Tally.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Cell_Tally::Cell_Tally(SP_Physics physics)
    : Base(physics, false)
    , d_geometry(b_physics->get_geometry())
    , d_state_idx(0)
{
    REQUIRE(d_geometry);

    // Add tally state metadata to the particle

    // First make a member memory manager
    auto member =
        std::make_shared<metaclass::Member_Manager_Cell_Tally_State>();

    // Now make the metadata in the particle
    const_cast<unsigned int &>(d_state_idx) =
        Particle_t::Metadata::new_member<State_t>("cell_tally_state", member);
    CHECK(Particle_t::Metadata::name(d_state_idx) == "cell_tally_state");

    // Set the tally name
    set_name("cell");

    // Seset tally
    reset();
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Add a list of cells to tally.
 */
void Cell_Tally::set_cells(const std::vector<int> &cells)
{
    // Make a new tally results
    Result tally;

    // Iterate through the cells and insert them into the map
    for (const auto &cell : cells)
    {
        tally.insert({cell, {0.0, 0.0}});
    }

    // Swap with the existing tally
    std::swap(tally, d_tally);
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*
 * \brief Accumulate first and second moments.
 */
void Cell_Tally::end_history(const Particle_t &p)
{
    // Get the particles private tally
    const State_t::History_Tally &hist =
        p.metadata().access<State_t>(d_state_idx).state();

    // Iterate through local tally results and add them to the permanent
    // results
    for (const auto &t : hist)
    {
        CHECK(d_tally.find(t.first) != d_tally.end());

        // Get the persistent tally result
        auto &r = d_tally[t.first];

        // Store the moments
	r.atomic_sum_into( t.second, t.second*t.second );
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Do post-processing on first and second moments.
 */
void Cell_Tally::finalize(double num_particles)
{
    REQUIRE(num_particles > 1);

    // Do a global reduction on moments
    std::vector<double> first(d_tally.size(),  0.0);
    std::vector<double> second(d_tally.size(), 0.0);

    int ctr = 0;
    for (const auto &t : d_tally)
    {
        // Get a reference to the moments
        auto &moments = t.second;

        // Copy the moments into the communciation buffer
        first[ctr]  = moments.first();
        second[ctr] = moments.second();

        // Update the counter
        ++ctr;
    }
    CHECK(ctr == d_tally.size());

    // Do global reductions on the moments
    profugus::global_sum(first.data(),  first.size());
    profugus::global_sum(second.data(), second.size());

    // Reset the counter
    ctr = 0;

    // Iterate through tally cells and build the variance and mean
    for (auto &t : d_tally)
    {
        CHECK(d_geometry->cell_volume(t.first) > 0.0);

        // Get the volume for the cell
        double inv_V = 1.0 / d_geometry->cell_volume(t.first);

        // Store 1/N
        double inv_N = 1.0 / static_cast<double>(num_particles);

        // Calculate means for this cell
        double avg_l  = first[ctr] * inv_N;
        double avg_l2 = second[ctr] * inv_N;

        // Get a reference to the moments
        auto &moments = t.second;

        // Store the sample mean
        moments.first() = avg_l * inv_V;

        // Calculate the variance
        double var = num_particles / (num_particles - 1) * inv_V * inv_V *
                     (avg_l2 - avg_l * avg_l);

        // Store the error of the sample mean
        moments.second() = std::sqrt(var * inv_N);

        // Update the counter
        ++ctr;
    }
    CHECK(ctr == d_tally.size());
}

//---------------------------------------------------------------------------//
/*
 * \brief Clear/re-initialize all tally values between solves.
 */
void Cell_Tally::reset()
{
    // Clear all current tally results (but keep existing cells in place)
    for (auto &t : d_tally)
    {
	t.second.first() = 0.0;
	t.second.second() = 0.0;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Output results.
 */
void Cell_Tally::output(const std::string &out_file)
{
#ifdef USE_HDF5

    profugus::Serial_HDF5_Writer writer;
    writer.open(out_file, profugus::HDF5_IO::APPEND);

    // Allocate space for the output
    std::vector<int>    cells(d_tally.size());
    std::vector<double> results(d_tally.size() * 2, 0.0);

    // Fill the output
    int n = 0, ctr = 0;
    for (const auto &t : d_tally)
    {
        // Store the tally cell
        cells[n] = t.first;
        ++n;

        // Store the mean/error
        const auto &moments = t.second;

        // Mean
        results[ctr] = moments.first();
        ++ctr;

        // Error
        results[ctr] = moments.second();
        ++ctr;
    }
    CHECK(ctr == results.size());

    // Make a decomposition for the results
    HDF5_IO::Decomp d;
    d.ndims  = 2;
    d.global = {cells.size(), 2};
    d.local  = {cells.size(), 2};
    d.offset = {0, 0};
    d.order  = HDF5_IO::ROW_MAJOR;

    // >>> WRITE THE OUTPUT

    writer.begin_group("cell_tally");

    // Write the cells
    writer.write("cell_ids", cells);

    // Write the results
    writer.create_incremental_dataspace<double>("moments", d);
    writer.write_incremental_data("moments", d, results.data());

    writer.end_group();

    writer.close();

#endif
}

//---------------------------------------------------------------------------//
/*
 * \brief Track particle and tally..
 */
void Cell_Tally::accumulate( const double step,
			     Particle_t &p ) const
{
    REQUIRE(Particle_t::Metadata::name(d_state_idx) == "cell_tally_state");

    // Get the cell index
    int cell = d_geometry->cell(p.geo_state());

    // Get the particle local history tally
    auto &hist = *p.metadata().access<State_t>(d_state_idx).data();

    // O(1) check to see if we need to tally it
    if (d_tally.find(cell) != d_tally.end())
    {
        CHECK(d_tally.count(cell) == 1);

        // Tally for the history
        hist[cell] += p.wt() * step;
    }
}

} // end namespace profugus

//---------------------------------------------------------------------------//
// end of MC/mc/Cell_Tally.cc
//---------------------------------------------------------------------------//
