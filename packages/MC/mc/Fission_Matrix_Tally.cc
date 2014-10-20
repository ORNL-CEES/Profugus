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

#include <sstream>
#include <algorithm>

#include "harness/Warnings.hh"
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
    : Base(physics, true)
    , d_data(std::make_shared<FM_Data>())
{
    REQUIRE(!db.is_null());
    REQUIRE(db->isParameter("num_cycles"));
    REQUIRE(db->isSublist("fission_matrix_db"));
    REQUIRE(fm_mesh);
    REQUIRE(physics->get_geometry());
    REQUIRE(d_data);

    // set the name
    Base::set_name("fission_matrix");

    // make the data stored for the common tallies
    d_data->d_geometry      = physics->get_geometry();
    d_data->d_fm_mesh       = fm_mesh;
    d_data->d_cycle_ctr     = 0;
    d_data->d_tally_started = false;

    // size the denominator of the fission matrix tally
    d_data->d_denominator.resize(d_data->d_fm_mesh->num_cells());
    std::fill(d_data->d_denominator.begin(), d_data->d_denominator.end(), 0.0);

    // get the fission matrix sublist
    ParameterList_t &opt = db->sublist("fission_matrix_db");

    // get the number of cycles in the problem
    int num_cycles = db->get<int>("num_cycles");

    // set the starting cycle
    d_data->d_cycle_start = opt.get<int>("start_cycle", 0);

    // set the cycle output (by default, we don't output during transport by
    // setting to 1-past the last cycle)
    d_data->d_cycle_out = opt.get<int>("output_cycle", num_cycles);

    // make the source and pathlength tallies
    b_pl_tally  = std::make_shared<PL_Tally>(physics, d_data);
    b_src_tally = std::make_shared<Src_Tally>(physics, d_data);

    // get the output problem name for the hdf5 diagnostic file
    if (d_data->d_cycle_out < num_cycles)
    {
#ifdef USE_HDF5
        CHECK(db->isParameter("problem_name"));
        d_filename = db->get<std::string>("problem_name") + "_fm.h5";

        // make the initial file
        d_writer.open(d_filename);
        d_writer.close();
#else
        ADD_WARNING("HDF5 not available in this build, turning fission "
                    << "matrix output off.");
#endif
    }

    // make the sparse matrix using the optimal number of initial buckets and
    // with the correct block size in the hasher
    Fission_Matrix_Processor::Idx_Hash h(d_data->d_fm_mesh->num_cells());
    Sparse_Matrix m(d_data->d_numerator.bucket_count(), h);
    std::swap(m, d_data->d_numerator);

    VALIDATE(d_data->d_cycle_out >= d_data->d_cycle_start,
             "Fission matrix tallying starting on cycle "
             << d_data->d_cycle_start << ", but output requested on "
             << d_data->d_cycle_out);
    ENSURE(d_data->d_geometry);
    ENSURE(d_data->d_fm_mesh);
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Manually build the global fission matrix using the current state of
 * the tally.
 *
 * The client can access (and reset) the global fission matrix through the
 * processor by accessing the processor() function.
 */
void Fission_Matrix_Tally::build_matrix()
{
    // return if we haven't started tallying, as we assume that all domains
    // will have entries in the fission matrix (for replicated) in order for
    // this to work
    if (!d_data->d_tally_started)
        return;

    // use the processor to build the global fission matrix
    d_processor.build_matrix(d_data->d_numerator, d_data->d_denominator);
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief End cycle processing of fission matrix.
 */
void Fission_Matrix_Tally::end_cycle(double num_particles)
{
#ifdef USE_HDF5

    // build the sparse-stored, ordered fission matrix if we are past the
    // output cycle
    if (d_data->d_cycle_ctr >= d_data->d_cycle_out)
    {
        // use the processor to build the global fission matrix
        d_processor.build_matrix(d_data->d_numerator, d_data->d_denominator);

        // open the file - writing is only on proc 0
        d_writer.open(d_filename, HDF5_IO::APPEND, 0);

        // make the cycle group
        std::ostringstream m;
        m << "cycle_" << d_data->d_cycle_ctr;
        d_writer.begin_group(m.str());

        // write the full matrix dimensions
        d_writer.write("size",
                       static_cast<int>(d_data->d_fm_mesh->num_cells()));

        // get the indices of the matrix and the matrix elements
        const auto &indices = d_processor.graph();
        const auto &matrix  = d_processor.matrix();
        CHECK(indices.size() == matrix.size());

        // write the number of non-zero elements
        d_writer.write("non_zero", static_cast<int>(indices.size()));

        // get an HDF5 decomposition for the indices
        HDF5_IO::Decomp d(2);
        d.global[0] = indices.size();
        d.global[1] = 2;
        d.order     = HDF5_IO::ROW_MAJOR;

        // make space for the global data
        d_writer.create_incremental_dataspace<int>("indices", d);

        // we will do a single write
        d.local[0] = d.global[0];
        d.local[1] = d.global[1];

        // write the indices
        d_writer.write_incremental_data("indices", d, &indices[0].first);

        // write the elements of the fission matrix
        d_writer.write("matrix", matrix);

        // clear memory in the processor
        d_processor.reset();

        // close out writing
        d_writer.end_group();
        d_writer.close();

        CHECK(indices.empty());
        CHECK(matrix.empty());
    }
#endif

    // update internal counter
    ++d_data->d_cycle_ctr;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reset the tallies.
 */
void Fission_Matrix_Tally::reset()
{
    // clear all the tallies
    std::fill(d_data->d_denominator.begin(), d_data->d_denominator.end(), 0.0);

    // reset tally starting
    d_data->d_tally_started = false;

    // reset the processor
    d_processor.reset();

    // make a new sparse mesh
    Sparse_Matrix n;
    Fission_Matrix_Processor::Idx_Hash h(d_data->d_fm_mesh->num_cells());
    Sparse_Matrix m(n.bucket_count(), h);
    std::swap(m, d_data->d_numerator);

    ENSURE(d_data->d_numerator.empty());
}

//---------------------------------------------------------------------------//
// SRC_TALLY DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Setup fission matrix tallying at birth.
 */
void Fission_Matrix_Tally::Src_Tally::birth(const Particle_t &p)
{
    // return if we haven't started tallying yet
    if (d_data->d_cycle_start > d_data->d_cycle_ctr)
        return;

    // set tally started flag
    d_data->d_tally_started = true;

    REQUIRE(p.metadata().name(d_data->d_birth_idx) == "fm_birth_cell");

    // get the particle's geometric state
    const auto &geo_state = p.geo_state();

    // initialize the particle in the fission-matrix mesh
    d_data->d_fm_mesh->initialize(d_data->d_geometry->position(geo_state),
                                  d_data->d_geometry->direction(geo_state),
                                  d_fm_state);

    // determine the birth cell in the fission matrix mesh
    auto mesh_idx = d_data->d_fm_mesh->cell(d_fm_state);
    CHECK(mesh_idx >= 0 && mesh_idx < d_data->d_fm_mesh->num_cells());

    // set it in the particle's metadata
    const_cast<Particle_t &>(p).metadata().access<int>(d_data->d_birth_idx)
        = mesh_idx;

    // tally the birth weight
    d_data->d_denominator[mesh_idx] += p.wt();
}

//---------------------------------------------------------------------------//
// PL_TALLY DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Accumulate fission matrix pathlength tally.
 */
void Fission_Matrix_Tally::PL_Tally::accumulate(double            step,
                                                const Particle_t &p)
{
    using geometry::INSIDE;

    // return if we haven't started tallying yet
    if (d_data->d_cycle_start > d_data->d_cycle_ctr)
        return;

    REQUIRE(p.metadata().name(d_data->d_birth_idx) == "fm_birth_cell");

    // get the particle's geometric state
    const auto &geo_state = p.geo_state();

    // initialize the particle in the fission-matrix mesh
    d_data->d_fm_mesh->initialize(d_data->d_geometry->position(geo_state),
                                  d_data->d_geometry->direction(geo_state),
                                  d_fm_state);

    // remaining distance before finishing
    double remaining = step;

    // current fission matrix mesh cell (i)
    int i = 0;

    // cell particle was born in (j)
    int j = p.metadata().access<int>(d_data->d_birth_idx);

    // distance of current step through the fission matrix
    double d = 0.0;

    // get weighted contribution to the fission matrix tally that is constant
    // across the step
    double keff = p.wt() * b_physics->total(physics::NU_FISSION, p);

    // track through the fission matrix and accumulate the fission matrix
    // elements
    while (remaining > 0.0 &&
           d_data->d_fm_mesh->boundary_state(d_fm_state) == INSIDE)
    {
        // calculate next step
        d = d_data->d_fm_mesh->distance_to_boundary(d_fm_state);
        if (d > remaining)
            d = remaining;

        // get the current cell
        i = d_data->d_fm_mesh->cell(d_fm_state);
        CHECK(i >= 0 && i < d_data->d_fm_mesh->num_cells());

        // tally the fission matrix contribution to the (i,j) element
        d_data->d_numerator[Idx(i, j)] += d * keff;

        // subtract this step from the reaming distance
        remaining -= d;

        // move the particle to the next cell boundary
        d_data->d_fm_mesh->move_to_surface(d_fm_state);
    }
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Tally.cc
//---------------------------------------------------------------------------//
