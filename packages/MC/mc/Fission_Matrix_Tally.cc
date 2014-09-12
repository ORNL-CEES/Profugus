//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Matrix_Tally.cc
 * \author Thomas M. Evans
 * \date   Tue Jul 22 15:09:01 2014
 * \brief  Fission_Matrix_Tally member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fission_Matrix_Tally.hh"

#include <algorithm>

#include "Particle.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Fission_Matrix_Tally::Fission_Matrix_Tally(RCP_Std_DB       db,
                                           SP_Physics       physics,
                                           SP_Mesh_Geometry fm_mesh)
    : Base(physics)
    , d_geometry(physics->get_geometry())
    , d_fm_mesh(fm_mesh)
    , d_numerator()
    , d_denominator(d_fm_mesh->num_cells(), 0.0)
    , d_birth_idx(Particle::Metadata::new_pod_member<int>("fm_birth_cell"))
    , d_cycle_ctr(0)
{
    REQUIRE(!db.is_null());
    REQUIRE(db->isParameter("num_cycles"));
    REQUIRE(db->isSublist("fission_matrix_db"));

    // get the fission matrix sublist
    ParameterList_t &opt = db->sublist("fission_matrix_db");

    // get the number of cycles in the problem
    int num_cycles = db->get<int>("num_cycles");

    // set the cycle output
    d_cycle_out = opt.get<int>("output_cycle", num_cycles - 1);

    // make the sparse matrix using the optimal number of initial buckets and
    // with the correct block size in the hasher
    Fission_Matrix_Processor::Idx_Hash h(d_fm_mesh->num_cells());
    Sparse_Matrix m(d_numerator.bucket_count(), h);
    std::swap(m, d_numerator);

    ENSURE(d_cycle_out < num_cycles);
    ENSURE(d_geometry);
    ENSURE(d_fm_mesh);
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Setup fission matrix tallying at birth.
 */
void Fission_Matrix_Tally::birth(const Particle_t &p)
{
    REQUIRE(p.metadata().name(d_birth_idx) == "fm_birth_cell");

    // get the particle's geometric state
    const auto &geo_state = p.geo_state();

    // initialize the particle in the fission-matrix mesh
    d_fm_mesh->initialize(d_geometry->position(geo_state),
                          d_geometry->direction(geo_state), d_fm_state);

    // determine the birth cell in the fission matrix mesh
    auto mesh_idx = d_fm_mesh->cell(d_fm_state);
    CHECK(mesh_idx >= 0 && mesh_idx < d_fm_mesh->num_cells());

    // set it in the particle's metadata
    const_cast<Particle_t &>(p).metadata().access<int>(d_birth_idx) = mesh_idx;

    // tally the birth weight
    d_denominator[mesh_idx] += p.wt();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Accumulate fission matrix tally.
 */
void Fission_Matrix_Tally::accumulate(double            step,
                                      const Particle_t &p)
{
    using geometry::INSIDE;

    REQUIRE(p.metadata().name(d_birth_idx) == "fm_birth_cell");

    // get the particle's geometric state
    const auto &geo_state = p.geo_state();

    // initialize the particle in the fission-matrix mesh
    d_fm_mesh->initialize(d_geometry->position(geo_state),
                          d_geometry->direction(geo_state), d_fm_state);

    // remaining distance before finishing
    double remaining = step;

    // current fission matrix mesh cell (i)
    int i = 0;

    // cell particle was born in (j)
    int j = p.metadata().access<int>(d_birth_idx);

    // distance of current step through the fission matrix
    double d = 0.0;

    // get weighted contribution to the fission matrix tally that is constant
    // across the step
    double keff = p.wt() * b_physics->total(physics::NU_FISSION, p);

    // track through the fission matrix and accumulate the fission matrix
    // elements
    while (remaining > 0.0 && d_fm_mesh->boundary_state(d_fm_state) == INSIDE)
    {
        // calculate next step
        d = d_fm_mesh->distance_to_boundary(d_fm_state);
        if (d > remaining)
            d = remaining;

        // get the current cell
        i = d_fm_mesh->cell(d_fm_state);
        CHECK(i >= 0 && i < d_fm_mesh->num_cells());

        // tally the fission matrix contribution to the (i,j) element
        d_numerator[Idx(i, j)] += d * keff;

        // subtract this step from the reaming distance
        remaining -= d;

        // move the particle to the next cell boundary
        d_fm_mesh->move_to_surface(d_fm_state);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief End cycle processing of fission matrix.
 */
void Fission_Matrix_Tally::end_cycle(double num_particles)
{
    // update internal counter
    ++d_cycle_ctr;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reset the tallies.
 */
void Fission_Matrix_Tally::reset()
{
    // clear all the tallies
    std::fill(d_denominator.begin(), d_denominator.end(), 0.0);

    // make a new sparse mesh
    Sparse_Matrix n;
    Fission_Matrix_Processor::Idx_Hash h(d_fm_mesh->num_cells());
    Sparse_Matrix m(n.bucket_count(), h);
    std::swap(m, d_numerator);

    ENSURE(d_numerator.empty());
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Tally.cc
//---------------------------------------------------------------------------//
