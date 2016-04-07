//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/General_Source.t.hh
 * \author Steven Hamilton
 * \date   Mon Apr 04 20:38:12 2016
 * \brief  General_Source template member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_General_Source_t_hh
#define mc_General_Source_t_hh

#include "comm/global.hh"
#include "General_Source.hh"
#include "Sampler.hh"
#include "Global_RNG.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param db
 * \param geometry
 * \param physics
 * \param rng_control
 */
template <class Geometry>
General_Source<Geometry>::General_Source(RCP_ParameterList  db,
                                         SP_Geometry        geometry,
                                         SP_Physics         physics,
                                         SP_RNG_Control     rng_control)
    : Base(geometry, physics, rng_control)
    , d_np_total(0)
    , d_np_domain(0)
    , d_wt(1.0)
    , d_np_left(0)
    , d_np_run(0)
{
    REQUIRE( Teuchos::nonnull(db) );
    REQUIRE( b_geometry );
    REQUIRE( b_physics );
    REQUIRE( b_rng_control );
    REQUIRE( b_nodes > 0 );

    // Resize CDFs
    int num_cells  = b_geometry->num_cells();
    int num_groups = b_physics->num_groups();
    d_cell_cdf.resize(num_cells,0.0);
    d_erg_cdfs.resize(num_cells,std::vector<double>(num_groups,0.0));

    int np_requested = static_cast<def::size_type>(db->get("Np",1000));
    d_np_domain = np_requested / b_nodes;

    // Recompute total (may change slightly from requested)
    d_np_total = d_np_domain * b_nodes;
    d_wt = static_cast<double>(np_requested) / static_cast<double>(d_np_total);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct source.
 */
template <class Geometry>
void General_Source<Geometry>::build_source(Teuchos::TwoDArray<double> &source)
{
    REQUIRE( source.getNumRows() == b_geometry->num_cells() );
    REQUIRE( source.getNumCols() == b_physics->num_groups() );

    double cell_sum = 0.0;
    for( int cell = 0; cell < b_geometry->num_cells(); ++cell )
    {
        const auto &cell_source = source[cell];
        CHECK( cell_source.size() == b_physics->num_groups() );

        // Compute total source for cell
        double this_cell = std::accumulate(cell_source.begin(),
                                           cell_source.end(),
                                           0.0);
        cell_sum += this_cell;
        d_cell_cdf[cell] = cell_sum;

        // Compute energy CDF for this cell
        auto &this_erg_cdf = d_erg_cdfs[cell];
        double erg_sum = 0.0;
        for( int g = 0; g < b_physics->num_groups(); ++g )
        {
            erg_sum += cell_source[g];
            this_erg_cdf[g] = erg_sum;
        }

        // Normalize energy CDF
        double norm_factor = 1.0 / this_erg_cdf.back();
        for( auto &val : this_erg_cdf )
            val *= norm_factor;
    }

    // Normalize cell CDF
    double norm_factor = 1.0 / d_cell_cdf.back();
    for( auto &val : d_cell_cdf )
        val *= norm_factor;

    d_np_left = d_np_domain;
    d_np_run  = 0;

    // Build RNG
    Base::make_RNG();

    profugus::global_barrier();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Generate particle from source
 */
template <class Geometry>
auto General_Source<Geometry>::get_particle() -> SP_Particle
{
    REQUIRE( d_np_left > 0 );
    REQUIRE( d_wt > 0.0 );
    REQUIRE( profugus::Global_RNG::d_rng.assigned() );

    using def::I;
    using def::J;
    using def::K;

    SP_Particle p;

    // Return null particle if no histories left
    if( d_np_left == 0 )
        return p;

    // Make particle
    p = std::make_shared<Particle_t>();
    p->set_rng(profugus::Global_RNG::d_rng);
    auto rng = p->rng();

    p->set_wt(d_wt);

    // Sample angle isotropically
    Space_Vector omega;
    Base::sample_angle(omega, rng);

    // Determine cell
    auto cell = sampler::sample_discrete_CDF(b_geometry->num_cells(),
                                             &d_cell_cdf[0],
                                             rng.ran());

    // Sample from cell's CDF to determine group
    auto g = sampler::sample_discrete_CDF(b_physics->num_groups(),
                                          &d_erg_cdfs[cell][0],
                                          rng.ran());
    p->set_group(g);

    // Now determine spatial location within cell
    Space_Vector r;
    auto box = b_geometry->get_cell_extents(cell);
    const auto &lower = box.lower();
    const auto &upper = box.upper();
    bool found = false;
    for( int itr = 0; itr < 1000; ++itr )
    {
        double x = lower[I] + rng.ran() * (upper[I] - lower[I]);
        double y = lower[J] + rng.ran() * (upper[J] - lower[J]);
        double z = lower[K] + rng.ran() * (upper[K] - lower[K]);

        Space_Vector r = {x, y, z};

        // Initialize particles geo state
        b_geometry->initialize(r, omega, p->geo_state());

        if( cell == b_geometry->cell( p->geo_state() ) )
        {
            found = true;
            break;
        }
    }
    ENSURE( found );

    p->live();

    // Update counters
    d_np_left--;
    d_np_run++;

    return p;
}

} // end namespace profugus

#endif // mc_General_Source_t_hh

//---------------------------------------------------------------------------//
//                 end of General_Source.t.hh
//---------------------------------------------------------------------------//
