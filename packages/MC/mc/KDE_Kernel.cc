//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/kde/KDE_Kernel.cc
 * \author Gregory Davidson
 * \date   Mon Feb 16 14:21:15 2015
 * \brief  KDE_Kernel class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "KDE_Kernel.hh"

#include "harness/DBC.hh"
#include "comm/global.hh"
#include "utils/Container_Functions.hh"
#include "Sampler.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTION
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
KDE_Kernel::KDE_Kernel(SP_Geometry geometry,
                       SP_Physics  physics,
                       double      coefficient,
                       double      exponent)
    : d_geometry(geometry)
    , d_physics(physics)
    , d_coefficient(coefficient)
    , d_exponent(exponent)
    , d_num_sampled(0)
    , d_num_accepted(0)
{
    REQUIRE(d_geometry);
    REQUIRE(d_physics);
    REQUIRE(d_coefficient > 0.0);
    REQUIRE(d_exponent > -1.0 && d_exponent < 0.0);

    // Setup the bandwidth map (start with a bandwidth of all zero)
    cell_type num_cells = geometry->num_cells();
    for (cell_type cellid = 0; cellid < num_cells; ++cellid)
    {
        d_bndwidth_map[cellid] = 0.0;
    }
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set the bandwidth for a cell.
 */
void KDE_Kernel::calc_bandwidths(
    const Physics::Fission_Site_Container &fis_sites)
{
    typedef geometry::cell_type  cell_type;

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
        cell_type cell = d_geometry->cell(fs);

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

        // Calculate the bandwidth
        d_bndwidth_map[cell] = d_coefficient*std::sqrt(variance)*
                               std::pow(N, d_exponent);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the bandwidth for a given cell.
 */
double KDE_Kernel::bandwidth(cell_type cellid) const
{
    REQUIRE(d_bndwidth_map.count(cellid) == 1);

    return d_bndwidth_map.find(cellid)->second;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the bandwidths for all cells.
 */
std::vector<double>
KDE_Kernel::get_bandwidths() const
{
    // Create the vector to hold the bandwidths
    std::vector<double> bandwidths(d_bndwidth_map.size());

    // Loop over the bandwidth map and fill the bandwidth vector
    std::vector<double>::iterator vec_iter = bandwidths.begin();
    for (const auto &map_elem : d_bndwidth_map)
    {
        CHECK(vec_iter != bandwidths.end());

        // Fill element and update iterator
        *vec_iter++ = map_elem.first;
    }

    return bandwidths;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a new position.
 *
 * If the new position is outside the fissionable region, it is rejected.
 */
KDE_Kernel::Space_Vector
KDE_Kernel::sample_position(const Space_Vector &orig_position,
                            RNG                &rng) const
{
    REQUIRE(d_physics);
    REQUIRE(d_geometry);
    REQUIRE(rng.assigned());

    // Get the cell at the position
    cell_type cellid = d_geometry->cell(orig_position);

    size_type failures = 0;
    do
    {
        // Sample Epanechnikov kernel
        double epsilon = sampler::sample_epan(rng);

        // Get the bandwidth
        CHECK(d_bndwidth_map.count(cellid) == 1);
        double bandwidth = d_bndwidth_map.find(cellid)->second;
        CHECK(bandwidth >= 0.0);

        // Create a new position
        Space_Vector new_pos(orig_position[def::X],
                             orig_position[def::Y],
                             orig_position[def::Z] + epsilon*bandwidth/2.0);

        // Get matid from sampled point (may raise error if outside
        // geometry)
        if (d_physics->is_fissionable(d_geometry->matid(new_pos)))
        {
            // Accept: sampled point is fissionable
            d_num_sampled += failures + 1;
            ++d_num_accepted;
            return new_pos;
        }

        // Increment failure counter.
        ++failures;
    } while (failures != 1000);

    // No luck
    d_num_sampled += failures;
    std::cout << "1000 consecutive nonfissionable rejections in KDE."
              << std::endl;
    return Space_Vector(0,0,0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the fraction of samples accepted inside the kernel.
 */
double KDE_Kernel::acceptance_fraction() const
{
    REQUIRE(d_num_sampled > 0);

    return (static_cast<double>(d_num_accepted) /
            static_cast<double>(d_num_sampled));
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Extract the local z-locations of the Space_Vector and broadcast to
 *        all cores.
 */
std::vector<KDE_Kernel::Space_Vector>
KDE_Kernel::communicate_sites(
    const Physics::Fission_Site_Container &fis_sites) const
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
                       [](const Physics::Fission_Site &fs) { return fs.r; });

        // Communicate the vector
        profugus::global_sum(&global_locs[0][def::X], 3*global_locs.size());
    }

    return global_locs;
}

//---------------------------------------------------------------------------//
} // end namespace shift

//---------------------------------------------------------------------------//
// end of MC/mc/kde/KDE_Kernel.cc
//---------------------------------------------------------------------------//
