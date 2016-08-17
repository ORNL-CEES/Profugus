//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Cell_Tally.t.hh
 * \author Thomas M. Evans
 * \date   Fri Aug 21 15:29:14 2015
 * \brief  Cell_Tally class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Cell_Tally_t_hh
#define MC_mc_Cell_Tally_t_hh

#include <cmath>
#include <algorithm>

#include "Utils/comm/global.hh"
#include "Cell_Tally.hh"
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
Cell_Tally<Geometry>::Cell_Tally(RCP_Std_DB db, SP_Physics physics)
    : Base(physics, false)
    , d_geometry(b_physics->get_geometry())
    , d_db(db)
{
    REQUIRE(d_geometry);

    // set the tally name
    set_name("cell");

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
void Cell_Tally<Geometry>::set_cells(const std::vector<int> &cells)
{
    // Make a new tally results
    Result tally;

    // Iterate through the cells and insert them into the map
    for (const auto &cell : cells)
    {
        VALIDATE(cell < d_geometry->num_cells(),
                "Cell tally index exceeds number of cells in geometry.");
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
template <class Geometry>
void Cell_Tally<Geometry>::end_history()
{
    // Iterate through local tally results and add them to the permanent
    // results
    for (const auto &t : d_hist)
    {
        CHECK(d_tally.find(t.first) != d_tally.end());

        // Get the persistent tally result
        auto &r = d_tally[t.first];

        // Store the moments
        r.first  += t.second;
        r.second += t.second * t.second;
    }

    // Clear the local tally
    clear_local();
}

//---------------------------------------------------------------------------//
/*
 * \brief Do post-processing on first and second moments.
 */
template <class Geometry>
void Cell_Tally<Geometry>::finalize(double num_particles)
{
    REQUIRE(num_particles > 1);

    // Do a global reduction on moments
    std::vector<int>    cells(d_tally.size(),  0.0);
    std::vector<double> first(d_tally.size(),  0.0);
    std::vector<double> second(d_tally.size(), 0.0);

    int ctr = 0;
    for (const auto &t : d_tally)
    {
        // Get a reference to the moments
        auto &moments = t.second;

        // Copy the moments into the communciation buffer
        cells[ctr]  = t.first;
        first[ctr]  = moments.first;
        second[ctr] = moments.second;

        // Update the counter
        ++ctr;
    }
    CHECK(ctr == d_tally.size());

    // Do global reductions on the moments
    profugus::global_sum(first.data(),  first.size());
    profugus::global_sum(second.data(), second.size());

    // Reset the counter
    ctr = 0;

    const auto &volumes = d_geometry->cell_volumes();

    // Iterate through tally cells and build the variance and mean
    for (auto &t : d_tally)
    {
        CHECK(volumes[t.first] > 0.0);

        // Get the volume for the cell
        double inv_V = 1.0 / volumes[t.first];

        // Store 1/N
        double inv_N = 1.0 / static_cast<double>(num_particles);

        // Calculate means for this cell
        double avg_l  = first[ctr] * inv_N;
        double avg_l2 = second[ctr] * inv_N;

        // Get a reference to the moments
        auto &moments = t.second;

        // Store the sample mean
        moments.first = avg_l * inv_V;

        // Calculate the variance
        double var = num_particles / (num_particles - 1) * inv_V * inv_V *
                     (avg_l2 - avg_l * avg_l);

        // Store the error of the sample mean
        moments.second = std::sqrt(var * inv_N);

        // Store values back into vectors for HDF5 writing
        first[ctr]  = moments.first;
        second[ctr] = moments.second;

        // Update the counter
        ++ctr;
    }
    CHECK(ctr == d_tally.size());

    // Determine permutation vector for sorted cells
    std::vector<std::pair<int,int>> sort_vec;
    for (int i = 0; i < cells.size(); ++i)
        sort_vec.push_back({cells[i],i});

    std::sort(sort_vec.begin(), sort_vec.end(),
        [](const std::pair<int,int> &lhs, const std::pair<int,int> &rhs)
        { return lhs.first < rhs.first; } );

    std::vector<int> cell_map(cells.size());
    for (int i = 0; i < sort_vec.size(); ++i)
        cell_map[i] = sort_vec[i].second;

    // Reorder vectors
    {
        std::vector<int>    tmp_cells  = cells;
        std::vector<double> tmp_first  = first;
        std::vector<double> tmp_second = second;
        for (int i = 0; i < cell_map.size(); ++i)
        {
            int ind = cell_map[i];
            cells[i]  = tmp_cells[ind];
            first[i]  = tmp_first[ind];
            second[i] = tmp_second[ind];
        }
    }

#ifdef USE_HDF5
    REQUIRE( d_db->isType<std::string>("problem_name") );
    std::string filename = d_db->get<std::string>("problem_name") + "_flux.h5";

    Serial_HDF5_Writer writer;
    writer.open(filename);
    writer.write("cells",cells);
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
void Cell_Tally<Geometry>::reset()
{
    // Clear the local tally
    clear_local();

    // Clear all current tally results (but keep existing cells in place)
    for (auto &t : d_tally)
    {
        t.second.first  = 0.0;
        t.second.second = 0.0;
    }

    ENSURE(d_hist.empty());
}

//---------------------------------------------------------------------------//
/*
 * \brief Track particle and tally..
 */
template <class Geometry>
void Cell_Tally<Geometry>::accumulate(double            step,
                                      const Particle_t &p)
{
    // Get the cell index
    int cell = d_geometry->cell(p.geo_state());

    // O(1) check to see if we need to tally it
    if (d_tally.find(cell) != d_tally.end())
    {
        CHECK(d_tally.count(cell) == 1);

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
void Cell_Tally<Geometry>::clear_local()
{
    // Clear the local tally
    History_Tally tally;
    std::swap(tally, d_hist);

    ENSURE(d_hist.empty());
}

} // end namespace profugus

#endif // MC_mc_Cell_Tally_t_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Cell_Tally.t.hh
//---------------------------------------------------------------------------//
