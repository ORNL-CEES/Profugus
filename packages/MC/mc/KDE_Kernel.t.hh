//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/KDE_Kernel.t.hh
 * \author Gregory Davidson
 * \date   Mon Feb 16 14:21:15 2015
 * \brief  KDE_Kernel class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_KDE_Kernel_t_hh
#define MC_mc_KDE_Kernel_t_hh

#include "KDE_Kernel.hh"

#include "harness/DBC.hh"
#include "comm/global.hh"
#include "utils/Container_Functions.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTION
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Geometry>
KDE_Kernel<Geometry>::KDE_Kernel(SP_Geometry geometry,
                                 SP_Physics  physics,
                                 double      coefficient,
                                 double      exponent)
    : b_geometry(geometry)
    , b_physics(physics)
    , b_coefficient(coefficient)
    , b_exponent(exponent)
    , b_num_sampled(0)
    , b_num_accepted(0)
{
    REQUIRE(b_geometry);
    REQUIRE(b_physics);
    REQUIRE(b_coefficient > 0.0);
    REQUIRE(b_exponent > -1.0 && b_exponent < 0.0);

    // Setup the bandwidth map (start with a bandwidth of all zero)
    cell_type num_cells = b_geometry->num_cells();
    for (cell_type cellid = 0; cellid < num_cells; ++cellid)
    {
        b_bndwidth_map[cellid] = 0.0;
    }
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set the bandwidth for a cell.
 */
template<class Geometry>
void KDE_Kernel<Geometry>::calc_bandwidths(
    const Fission_Site_Container &fis_sites)
{
    typedef geometry::cell_type  cell_type;

    // Reset all bandwidths to zero
    for (cell_type cellid = 0; cellid < b_geometry->num_cells(); ++cellid)
    {
        b_bndwidth_map[cellid] = 0.0;
    }

    // Broadcast all of the fission sites to all cores
    std::vector<Space_Vector> global_sites = this->communicate_sites(fis_sites);

    // Stores a map of the sum and sum-squared of the z-position for each
    // fission site in each cell as well as the number of fission sites in
    // each cell
    std::map<cell_type, double> cell_sums;
    std::map<cell_type, double> cell_sum_sq;
    std::map<cell_type, unsigned int> num_fs;
    for (const Space_Vector& fs : global_sites)
    {
        // Get the Z-position
        double z = fs[def::Z];

        // Get the cell
        cell_type cell = b_geometry->cell(fs);

        cell_sums[cell]   += z;
        cell_sum_sq[cell] += z*z;
        num_fs[cell]      += 1;
    }

    // Loop over the cells and calculate the bandwidths
    auto cell_sum_iter = cell_sums.cbegin();
    auto cell_sq_iter  = cell_sum_sq.cbegin();
    for (auto num_fs_iter = num_fs.cbegin(); num_fs_iter != num_fs.cend();
         ++num_fs_iter, ++cell_sum_iter, ++cell_sq_iter)
    {
        CHECK(cell_sum_iter != cell_sums.cend());
        CHECK(cell_sq_iter  != cell_sum_sq.cend());
        CHECK(cell_sum_iter->first == cell_sq_iter->first);
        CHECK(cell_sum_iter->first == num_fs_iter->first);

        // Get the data
        cell_type cell = cell_sum_iter->first;
        double sum     = cell_sum_iter->second;
        double sum_sq  = cell_sq_iter->second;
        unsigned int N = num_fs_iter->second;

        // Calculate the variance
        double variance = sum_sq/N - (sum/N)*(sum/N);

        // Low particle counts can results in very small variance, which can
        // go to zero with roundoff.  Change to zero if this happens
        variance = std::max(variance, 0.0);

        // Calculate the bandwidth
        b_bndwidth_map[cell] = b_coefficient*std::sqrt(variance)*
                               std::pow(N, b_exponent);
        CHECK(b_bndwidth_map[cell] >= 0.0);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the bandwidth for a given cell.
 */
template<class Geometry>
double KDE_Kernel<Geometry>::bandwidth(cell_type cellid) const
{
    REQUIRE(b_bndwidth_map.count(cellid) == 1);

    return b_bndwidth_map.find(cellid)->second;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the bandwidths for all cells.
 */
template<class Geometry>
std::vector<double>
KDE_Kernel<Geometry>::get_bandwidths() const
{
    // Create the vector to hold the bandwidths
    std::vector<double> bandwidths(b_bndwidth_map.size());

    // Loop over the bandwidth map and fill the bandwidth vector
    std::vector<double>::iterator vec_iter = bandwidths.begin();
    for (const auto &map_elem : b_bndwidth_map)
    {
        CHECK(vec_iter != bandwidths.end());

        // Fill element and update iterator
        *vec_iter++ = map_elem.first;
    }

    return bandwidths;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the bandwidth for a given cell.
 */
template<class Geometry>
void KDE_Kernel<Geometry>::set_bandwidth(geometry::cell_type cell,
                                         double              bandwidth)
{
    REQUIRE(b_bndwidth_map.find(cell) != b_bndwidth_map.end());

    b_bndwidth_map[cell] = bandwidth;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the fraction of samples accepted inside the kernel.
 */
template<class Geometry>
double KDE_Kernel<Geometry>::acceptance_fraction() const
{
    REQUIRE(b_num_sampled > 0);

    return (static_cast<double>(b_num_accepted) /
            static_cast<double>(b_num_sampled));
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Extract the local z-locations of the Space_Vector and broadcast to
 *        all cores.
 */
template<class Geometry>
std::vector<typename KDE_Kernel<Geometry>::Space_Vector>
KDE_Kernel<Geometry>::communicate_sites(
    const Fission_Site_Container &fis_sites) const
{
    // >>> COMMUNICATE SIZES
    // Create a vector to hold all of the local sizes
    std::vector<int> local_sizes(profugus::nodes(), 0);

    // Set the local size
    local_sizes[profugus::node()] = fis_sites.size();

    // Communicate all of the local sizes
    profugus::global_sum(local_sizes.data(), local_sizes.size());

    // Calculate the global size
    int global_size =
        std::accumulate(local_sizes.begin(), local_sizes.end(), 0);

    // >>> COMMUNICATE GLOBAL FISSION SITES
    // Create a vector to hold all of the fission sites
    std::vector<Space_Vector> global_locs(global_size,
                                          Space_Vector(0.0, 0.0, 0.0));
    if (global_size > 0)
    {
        // Find the copy index
        int begin_index =
            std::accumulate(local_sizes.begin(),
                            local_sizes.begin() + profugus::node(), 0);
        CHECK(begin_index + fis_sites.size() <= global_locs.size());

        // Fill the local portion with the z locations
        std::transform(fis_sites.begin(), fis_sites.end(),
                       global_locs.begin() + begin_index,
                       [](const typename Physics_t::Fission_Site &fs)
                       { return fs.r; });

        // Communicate the vector
        profugus::global_sum(&global_locs[0][def::X], 3*global_locs.size());
    }

    return global_locs;
}

//---------------------------------------------------------------------------//
} // end namespace profugus

#endif // MC_mc_KDE_Kernel_t_hh

//---------------------------------------------------------------------------//
// end of MC/mc/kde/KDE_Kernel.t.hh
//---------------------------------------------------------------------------//
