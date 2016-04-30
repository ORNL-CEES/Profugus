//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Current_Tally.t.hh
 * \author Steven Hamilton
 * \date   Thu Apr 28 20:19:42 2016
 * \brief  Current_Tally template member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Current_Tally_t_hh
#define mc_Current_Tally_t_hh

#include <algorithm>

#include "Current_Tally.hh"

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"
#include "utils/Vector_Functions.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Current_Tally<Geometry>::Current_Tally(RCP_ParameterList db,
                                       SP_Physics physics,
                                       const std::vector<double> &x_edges,
                                       const std::vector<double> &y_edges,
                                       const std::vector<double> &z_edges)
    : Base(physics, false)
    , d_x_edges(x_edges)
    , d_y_edges(y_edges)
    , d_z_edges(z_edges)
{
    REQUIRE(x_edges.size() > 1);
    REQUIRE(y_edges.size() > 1);
    REQUIRE(z_edges.size() > 1);

    int num_x_edges = x_edges.size();
    int num_y_edges = y_edges.size();
    int num_z_edges = z_edges.size();

    int num_x_faces = num_x_edges * (num_y_edges - 1) * (num_z_edges - 1);
    int num_y_faces = (num_x_edges - 1) * num_y_edges * (num_z_edges - 1);
    int num_z_faces = (num_x_edges - 1) * (num_y_edges - 1) * num_z_edges;

    d_x_current.resize(num_x_faces);
    d_y_current.resize(num_y_faces);
    d_z_current.resize(num_z_faces);
    d_x_std_dev.resize(num_x_faces);
    d_y_std_dev.resize(num_y_faces);
    d_z_std_dev.resize(num_z_faces);
    d_x_hist.resize(num_x_faces);
    d_y_hist.resize(num_y_faces);
    d_z_hist.resize(num_z_faces);

    // Compute surface areas
    d_x_areas.resize(num_x_faces,0.0);
    d_y_areas.resize(num_y_faces,0.0);
    d_z_areas.resize(num_z_faces,0.0);

    // X surfaces
    for (int i = 0; i < num_x_edges; ++i)
    {
        for (int j = 0; j < num_y_edges-1; ++j)
        {
            double y_width = y_edges[j+1] - y_edges[j];
            CHECK(y_width > 0.0);
            for (int k = 0; k < num_z_edges-1; ++k)
            {
                double z_width = z_edges[k+1] - z_edges[k];
                CHECK(z_width > 0.0);
                int ind = i + num_x_edges * (j + (num_y_edges - 1) * k);
                CHECK(ind < d_x_areas.size());
                d_x_areas[ind] = y_width * z_width;
            }
        }
    }

    // Y surfaces
    for (int i = 0; i < num_x_edges-1; ++i)
    {
        double x_width = x_edges[i+1] - x_edges[i];
        CHECK(x_width > 0.0);
        for (int j = 0; j < num_y_edges; ++j)
        {
            for (int k = 0; k < num_z_edges-1; ++k)
            {
                double z_width = z_edges[k+1] - z_edges[k];
                CHECK(z_width > 0.0);
                int ind = i + (num_x_edges - 1) * (j + num_y_edges * k);
                CHECK(ind < d_y_areas.size());
                d_y_areas[ind] = x_width * z_width;
            }
        }
    }

    // Z surfaces
    for (int i = 0; i < num_x_edges-1; ++i)
    {
        double x_width = x_edges[i+1] - x_edges[i];
        CHECK(x_width > 0.0);
        for (int j = 0; j < num_y_edges-1; ++j)
        {
            double y_width = y_edges[j+1] - y_edges[j];
            CHECK(y_width > 0.0);
            for (int k = 0; k < num_z_edges; ++k)
            {
                int ind = i + (num_x_edges - 1) * (j + (num_y_edges - 1) * k);
                CHECK(ind < d_z_areas.size());
                d_z_areas[ind] = x_width * y_width;
            }
        }
    }

    // set the tally name
    this->set_name("current");

    // reset tally
    reset();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Accumulate tally on surface
 */
template <class Geometry>
void Current_Tally<Geometry>::tally_surface(const Particle_t &particle)
{
    using def::I;
    using def::J;
    using def::K;
    using profugus::soft_equiv;

    constexpr double tol = 1e-12;

    const auto &xyz = particle.geo_state().d_r;

    // Locate particle on surface
    // Each check will override previous checks
    // We process these in z-y-x order so that if a particle is on an
    // edge or a corner, the precedence will be x-y-z
    int edge = -1;

    // Locate z region of particle
    int z_face = -1;
    auto itr = std::lower_bound(d_z_edges.begin(), d_z_edges.end(), xyz[K]);
    if (soft_equiv(*itr,xyz[K],tol))
    {
        edge = K;
    }
    else if (itr != d_z_edges.begin() && soft_equiv(*(--itr),xyz[K],tol))
    {
        edge = K;
    }
    z_face = itr - d_z_edges.begin();

    // Locate y region of particle
    int y_face = -1;
    itr = std::lower_bound(d_y_edges.begin(), d_y_edges.end(), xyz[J]);
    if (soft_equiv(*itr,xyz[J],tol))
    {
        edge = J;
    }
    else if (itr != d_y_edges.begin() && soft_equiv(*(--itr),xyz[J],tol))
    {
        edge = J;
    }
    y_face = itr - d_y_edges.begin();

    // Locate x region of particle
    int x_face = -1;
    itr = std::lower_bound(d_x_edges.begin(), d_x_edges.end(), xyz[I]);
    if (soft_equiv(*itr,xyz[I],tol))
    {
        edge = I;
    }
    else if (itr != d_x_edges.begin() && soft_equiv(*(--itr),xyz[I],tol))
    {
        edge = I;
    }
    x_face = itr - d_x_edges.begin();

    if (edge == I)
    {
        int ind = x_face + d_x_edges.size() *
            (y_face + (d_y_edges.size() - 1) * z_face);
        ENSURE (ind < d_x_hist.size());


        const auto &dir = particle.geo_state().d_dir;
        double dot = profugus::dot_product(dir, {1.0, 0.0, 0.0});

        d_x_hist[ind] += particle.wt() * dot;
    }
    else if (edge == J)
    {
        int ind = x_face + (d_x_edges.size() - 1) *
            (y_face + d_y_edges.size() * z_face);
        ENSURE (ind < d_y_hist.size());

        const auto &dir = particle.geo_state().d_dir;
        double dot = profugus::dot_product(dir, {0.0, 1.0, 0.0});

        d_y_hist[ind] += particle.wt() * dot;
    }
    else if (edge == K)
    {
        int ind = x_face + (d_x_edges.size() - 1) *
            (y_face + (d_y_edges.size() - 1) * z_face);
        ENSURE (ind < d_z_hist.size());

        const auto &dir = particle.geo_state().d_dir;
        double dot = profugus::dot_product(dir, {0.0, 0.0, 1.0});

        d_z_hist[ind] += particle.wt() * dot;
    }

}

//---------------------------------------------------------------------------//
/*!
 * \brief Finalize tally result
 */
template <class Geometry>
void Current_Tally<Geometry>::end_history()
{
    REQUIRE( d_x_hist.size() == d_x_current.size() );
    REQUIRE( d_x_hist.size() == d_x_std_dev.size() );
    REQUIRE( d_y_hist.size() == d_y_current.size() );
    REQUIRE( d_y_hist.size() == d_y_std_dev.size() );
    REQUIRE( d_z_hist.size() == d_z_current.size() );
    REQUIRE( d_z_hist.size() == d_z_std_dev.size() );

    // Add contribution to mean and variance
    for (int i = 0; i < d_x_hist.size(); ++i)
    {
        d_x_current[i] += d_x_hist[i];
        d_x_std_dev[i] += d_x_hist[i]*d_x_hist[i];
    }
    for (int j = 0; j < d_y_hist.size(); ++j)
    {
        d_y_current[j] += d_y_hist[j];
        d_y_std_dev[j] += d_y_hist[j]*d_y_hist[j];
    }
    for (int k = 0; k < d_z_hist.size(); ++k)
    {
        d_z_current[k] += d_z_hist[k];
        d_z_std_dev[k] += d_z_hist[k]*d_z_hist[k];
    }

    // Reset history vectors
    std::fill(d_x_hist.begin(), d_x_hist.end(), 0.0);
    std::fill(d_y_hist.begin(), d_y_hist.end(), 0.0);
    std::fill(d_z_hist.begin(), d_z_hist.end(), 0.0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Finalize tally result
 */
template <class Geometry>
void Current_Tally<Geometry>::finalize(double num_particles)
{
    REQUIRE( d_x_std_dev.size() == d_x_current.size() );
    REQUIRE( d_y_std_dev.size() == d_y_current.size() );
    REQUIRE( d_z_std_dev.size() == d_z_current.size() );

    // Do global reduction on first and second moments
    profugus::global_sum(d_x_current.data(), d_x_current.size());
    profugus::global_sum(d_y_current.data(), d_y_current.size());
    profugus::global_sum(d_z_current.data(), d_z_current.size());
    profugus::global_sum(d_x_std_dev.data(), d_x_std_dev.size());
    profugus::global_sum(d_y_std_dev.data(), d_y_std_dev.size());
    profugus::global_sum(d_z_std_dev.data(), d_z_std_dev.size());

    double inv_N = 1.0 / num_particles;

    // Finalize x currents and std deviations
    for (int i = 0; i < d_x_current.size(); ++i)
    {
        double inv_A = 1.0 / d_x_areas[i];

        double avg_l  = d_x_current[i] * inv_N;
        double avg_l2 = d_x_std_dev[i] * inv_N;

        d_x_current[i] = avg_l * inv_A;

        double var = num_particles / (num_particles - 1) * inv_A * inv_A *
            (avg_l2 - avg_l * avg_l);

        d_x_std_dev[i] = std::sqrt(var * inv_N);
    }

    // Finalize y currents and std deviations
    for (int i = 0; i < d_y_current.size(); ++i)
    {
        double inv_A = 1.0 / d_y_areas[i];

        double avg_l  = d_y_current[i] * inv_N;
        double avg_l2 = d_y_std_dev[i] * inv_N;

        d_y_current[i] = avg_l * inv_A;

        double var = num_particles / (num_particles - 1) * inv_A * inv_A *
            (avg_l2 - avg_l * avg_l);

        d_y_std_dev[i] = std::sqrt(var * inv_N);
    }

    // Finalize z currents and std deviations
    for (int i = 0; i < d_z_current.size(); ++i)
    {
        double inv_A = 1.0 / d_z_areas[i];

        double avg_l  = d_z_current[i] * inv_N;
        double avg_l2 = d_z_std_dev[i] * inv_N;

        d_z_current[i] = avg_l * inv_A;

        double var = num_particles / (num_particles - 1) * inv_A * inv_A *
            (avg_l2 - avg_l * avg_l);

        d_z_std_dev[i] = std::sqrt(var * inv_N);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reset tally
 */
template <class Geometry>
void Current_Tally<Geometry>::reset()
{
    std::fill(d_x_current.begin(), d_x_current.end(), 0.0);
    std::fill(d_y_current.begin(), d_y_current.end(), 0.0);
    std::fill(d_z_current.begin(), d_z_current.end(), 0.0);
    std::fill(d_x_std_dev.begin(), d_x_std_dev.end(), 0.0);
    std::fill(d_y_std_dev.begin(), d_y_std_dev.end(), 0.0);
    std::fill(d_z_std_dev.begin(), d_z_std_dev.end(), 0.0);
    std::fill(d_x_hist.begin(), d_x_hist.end(), 0.0);
    std::fill(d_y_hist.begin(), d_y_hist.end(), 0.0);
    std::fill(d_z_hist.begin(), d_z_hist.end(), 0.0);
}

} // end namespace profugus

#endif // mc_Current_Tally_t_hh

//---------------------------------------------------------------------------//
//                 end of Current_Tally.t.hh
//---------------------------------------------------------------------------//
