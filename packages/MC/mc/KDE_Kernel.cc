//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Shift/mc_sources/kde/KDE_Kernel.cc
 * \author Seth R Johnson
 * \date   Mon Feb 16 14:21:15 2015
 * \brief  KDE_Kernel class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "KDE_Kernel.hh"

#include "Nemesis/harness/DBC.hh"
#include "Nemesis/comm/global.hh"
#include "Nemesis/utils/Container_Functions.hh"

namespace shift
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
void Axial_Kernel::calc_bandwidths(const std::vector<Space_Vector> &fis_sites)
{
    typedef geometria::cell_type  cell_type;

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
        cell_type cell = b_geometry->find_cell(fs);

        cell_sums[cell]   += z;
        cell_sum_sq[cell] += z*z;
        num_fs[cell]      += 1;
    }

    // Loop over the cells and calculate the bandwidths
    std::map<cell_type, double> bandwidths;
    for (const auto &zip_elem : nemesis::zip(cell_sums, cell_sum_sq, num_fs))
    {
        cell_type cell = std::get<0>(std::get<0>(zip_elem));
        double sum     = std::get<1>(std::get<0>(zip_elem));
        double sum_sq  = std::get<1>(std::get<1>(zip_elem));
        unsigned int N = std::get<1>(std::get<2>(zip_elem));

        // Calculate the variance
        double variance = sum_sq/N - (sum/N)*(sum/N);

        // Calculate the bandwidth
        bandwidths[cell] = b_coefficient*std::sqrt(variance)*
                           std::pow(N, b_exponent);
    }

    // Calculate the bandwidth
    d_bandwidth_map.set(bandwidths);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the bandwidth for a given cell.
 */
double Axial_Kernel::bandwidth(cell_type cellid) const
{
    REQUIRE(d_bndwidth_map.count(cellid) == 1);

    return d_bndwidth_map[cellid];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a new position.
 *
 * If the new position is outside the fissionable region, it is rejected.
 */
Axial_Kernel::Space_Vector
Axial_Kernel::sample_position(const Space_Vector &orig_position,
                              RNG                &rng) const
{
    REQUIRE(b_physics);
    REQUIRE(b_geometry);
    REQUIRE(rng.assigned());

    // Get the cell at the position
    cell_type cellid = b_geometry->find_cell(orig_position);

    size_type failures = 0;
    do
    {
        // Sample Epanechnikov kernel
        double epsilon = mc::sampler::sample_epan(rng);

        // If we have a bandwidth here, get it.  Otherwise, use a bandwidth of
        // 0.0
        double bandwidth = 0.0;
        if (d_bandwidth_map.have_cellid(cellid))
        {
            bandwidth = d_bndwidth_map.bandwidth(cellid);
        }
        Check(bandwidth >= 0.0);

        // Create a new position
        Space_Vector new_pos(orig_position[def::X],
                             orig_position[def::Y],
                             orig_position[def::Z] + epsilon*bandwidth/2.0);

        try
        {
            // Get matid from sampled point (may raise error if outside
            // geometry)
            if (b_physics->is_fissionable(b_geometry->find_matid(new_pos)))
            {
                // Accept: sampled point is fissionable
                b_num_sampled += failures + 1;
                ++b_num_accepted;
                return new_pos;
            }
        }
        catch (const geometria::Geometry_Error& e)
        {
            // Outside the geometry!
            continue;
        }

        // Increment failure counter.
        ++failures;
    } while (failures != 1000);

    // No luck
    b_num_sampled += failures;
    throw geometria::Geometry_Error(
        "1000 consecutive nonfissionable rejections in KDE.");
    return Space_Vector(0,0,0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the fraction of samples accepted inside the kernel.
 */
double KDE_Kernel::acceptance_fraction() const
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
std::vector<KDE_Kernel::Space_Vector>
KDE_Kernel::communicate_sites(const std::vector<Space_Vector> &fis_sites) const
{
    // >>> COMMUNICATE SIZES
    // Create a vector to hold all of the local sizes
    std::vector<int> local_sizes(nemesis::nodes(), 0);

    // Set the local size
    local_sizes[nemesis::node()] = fis_sites.size();

    // Communicate all of the local sizes
    nemesis::global_sum(local_sizes.data(), local_sizes.size());

    // Calculate the global size
    int global_size = nemesis::accumulate(local_sizes);

    // >>> COMMUNICATE GLOBAL Z-LOCATIONS
    // Create a vector to hold all of the z-locations
    std::vector<Space_Vector> global_locs(global_size,
                                          Space_Vector(0.0, 0.0, 0.0));

    if (global_size > 0)
    {
        // Find the copy index
        int begin_index =
            std::accumulate(local_sizes.begin(),
                            local_sizes.begin() + nemesis::node(), 0);
        CHECK(begin_index + fis_sites.size() <= global_locs.size());

        // Fill the local portion with the z locations
        std::copy(fis_sites.begin(), fis_sites.end(),
                  global_locs.begin() + begin_index);

        // Communicate the vector
        nemesis::global_sum(&global_locs[0][def::X], 3*global_locs.size());
    }

    return global_locs;
}

//---------------------------------------------------------------------------//
} // end namespace shift

//---------------------------------------------------------------------------//
// end of Shift/mc_sources/kde/KDE_Kernel.cc
//---------------------------------------------------------------------------//
