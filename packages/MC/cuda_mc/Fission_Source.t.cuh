//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fission_Source.t.cuh
 * \author Stuart Slattery
 * \date   Mon May 05 14:22:46 2014
 * \brief  Fission_Source template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Fission_Source_t_cuh
#define cuda_mc_Fission_Source_t_cuh

#include <algorithm>
#include <numeric>

#include "Teuchos_Array.hpp"

#include "Fission_Source.hh"

#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Memory.cuh"
#include "cuda_utils/Hardware.hh"
#include "cuda_utils/Utility_Functions.hh"
#include "cuda_utils/Definitions.hh"
#include "comm/global.hh"
#include "mc/Global_RNG.hh"

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
// CUDA KERNELS
//---------------------------------------------------------------------------//
// Sample the mesh for the initial particle state.
template<class Geometry>
__global__ void sample_mesh_kernel( const Geometry* geometry,
				    const Physics<Geometry>* physics,
				    const Cartesian_Mesh* fission_mesh,
				    const int* fission_cells,
				    const int num_particle,
				    const double weight,
                                    const int active_cycle,
				    Particle_Vector<Geometry>* particles )
{
    using def::I; using def::J; using def::K;

    // Get the thread index.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int start_idx = particles->event_lower_bound( events::DEAD );

    if ( idx < num_particle )
    {
	// Get the particle index.
	int pidx = idx + start_idx;

	// sample the angle isotropically
	cuda_utils::Space_Vector omega;
	cuda_utils::utility::sample_angle( 
	    omega, particles->ran(pidx), particles->ran(pidx) );

	// get the logical indices of the birth cell
	int idim, jdim, kdim;
	fission_mesh->cardinal( fission_cells[idx], idim, jdim, kdim );

        // get the low/high bounds for the cell
        double xdims[] = {fission_mesh->edges(I)[idim],
                          fission_mesh->edges(I)[idim + 1]};
        double ydims[] = {fission_mesh->edges(J)[jdim],
                          fission_mesh->edges(J)[jdim + 1]};
        double zdims[] = {fission_mesh->edges(K)[kdim],
                          fission_mesh->edges(K)[kdim + 1]};

	// sample the mesh until a fission site is found
	int group = -1;
	int matid = -1;
	bool sampled = false;
	cuda_utils::Space_Vector r;
	while (!sampled)
	{
	    // sample a point in the geometry
	    r.x = (xdims[1] - xdims[0]) * particles->ran(pidx) + xdims[0];
	    r.y = (ydims[1] - ydims[0]) * particles->ran(pidx) + ydims[0];
	    r.z = (zdims[1] - zdims[0]) * particles->ran(pidx) + zdims[0];

	    // initialize the geometry state
	    geometry->initialize( r, omega, particles->geo_state(pidx) );

	    // get the material id
	    matid = geometry->matid( particles->geo_state(pidx) );

	    // try to initialize fission. exit the loop if successful.
	    physics->initialize_fission( matid,
					 particles->ran(pidx),
					 group,
					 sampled );
	}

	// set the matid
	particles->set_matid( pidx, matid );

	// set the group
	particles->set_group( pidx, group );

	// set the weight
	particles->set_wt( pidx, weight );

        // sample the initial distance to collision
        particles->set_dist_mfp( pidx, -std::log(particles->ran(pidx)) );

	// set the event.
	particles->set_event( pidx, events::TAKE_STEP );

        // set the batch.
        particles->set_batch( pidx, active_cycle );

	// make particle alive
	particles->live( pidx );
    }
}

//---------------------------------------------------------------------------//
// Sample the geometry for the initial particle state.
template<class Geometry>
__global__ void sample_geometry_kernel( const Geometry* geometry,
					const Physics<Geometry>* physics,
					const int num_particle,
					const double weight,
					const cuda_utils::Space_Vector width,
					const cuda_utils::Space_Vector lower,
                                        const int active_cycle,
					Particle_Vector<Geometry>* particles )
{
    // Get the thread index.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int start_idx = particles->event_lower_bound( events::DEAD );

    if ( idx < num_particle )
    {
	// Get the particle index.
	int pidx = idx + start_idx;

	// sample the angle isotropically
	cuda_utils::Space_Vector omega;
	cuda_utils::utility::sample_angle( 
	    omega, particles->ran(pidx), particles->ran(pidx) );

	// sample the geometry until a fission site is found
	int group = -1;
	int matid = -1;
	bool sampled = false;
	cuda_utils::Space_Vector r;
	while (!sampled)
	{
	    // sample a point in the geometry
	    r.x = width.x * particles->ran(pidx) + lower.x;
	    r.y = width.y * particles->ran(pidx) + lower.y;
	    r.z = width.z * particles->ran(pidx) + lower.z;

	    // initialize the geometry state
	    geometry->initialize( r, omega, particles->geo_state(pidx) );

	    // get the material id
	    matid = geometry->matid( particles->geo_state(pidx) );

	    // try to initialize fission. exit the loop if successful.
	    physics->initialize_fission( matid,
					 particles->ran(pidx),
					 group,
					 sampled );
	}

	// set the matid
	particles->set_matid( pidx, matid );

	// set the group
	particles->set_group( pidx, group );

	// set the weight
	particles->set_wt( pidx, weight );

        // sample the initial distance to collision
        particles->set_dist_mfp( pidx, -std::log(particles->ran(pidx)) );

	// set the event.
	particles->set_event( pidx, events::TAKE_STEP );

        // set the batch.
        particles->set_batch( pidx, active_cycle );

	// make particle alive
	particles->live( pidx );
    }
}

//---------------------------------------------------------------------------//
// Sample the fission sites for the initial 
template<class Geometry>
__global__ void sample_fission_sites_kernel(
    const Geometry* geometry,
    const Physics<Geometry>* physics,
    const int num_particle,
    const typename Physics<Geometry>::Fission_Site* fission_sites,
    const double weight,
    const int active_cycle,
    Particle_Vector<Geometry>* particles )
{
    // Get the thread index.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int start_idx = particles->event_lower_bound( events::DEAD );

    if ( idx < num_particle )
    {
	// Get the particle index.
	int pidx = idx + start_idx;

	// sample the angle isotropically
	cuda_utils::Space_Vector omega;
	cuda_utils::utility::sample_angle( 
	    omega, particles->ran(pidx), particles->ran(pidx) );

	// initialize the geometry state
	geometry->initialize( fission_sites[idx].r, 
			      omega, 
			      particles->geo_state(pidx) );

	// get the matid.
	int matid = geometry->matid( particles->geo_state(pidx) );

	// initialize fission
	bool sampled = false;
	int group = -1;
	physics->initialize_fission( matid, 
				     particles->ran(pidx),
				     group,
				     sampled );

	// set the matid
	particles->set_matid( pidx, matid );

	// set the group
	particles->set_group( pidx, group );

	// set the weight
	particles->set_wt( pidx, weight );

        // sample the initial distance to collision
        particles->set_dist_mfp( pidx, -std::log(particles->ran(pidx)) );

	// set the event.
	particles->set_event( pidx, events::TAKE_STEP );

        // set the batch.
        particles->set_batch( pidx, active_cycle );

	// make particle alive
	particles->live( pidx );
   }
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Fission_Source<Geometry>::Fission_Source(const RCP_Std_DB&     db,
                                         const SDP_Geometry&   geometry,
                                         const SDP_Physics&    physics)
    : d_geometry(geometry)
    , d_physics(physics)
    , d_fission_rebalance(std::make_shared<Fission_Rebalance_t>())
    , d_np_requested(0)
    , d_np_total(0)
    , d_np_domain(0)
    , d_wt(0.0)
    , d_np_left(0)
    , d_np_run(0)
    , d_active_cycle_counter(0)
{
    using def::I; using def::J; using def::K;

    REQUIRE(d_geometry);
    REQUIRE(d_physics);

    // Boundaries in -X, +X, -Y, +Y, -Z, +Z
    Teuchos::Array<double> extents(6, 0.);

    // Assign large extents that will be trimmed to the geometry by default
    extents[0] = extents[2] = extents[4] = -std::numeric_limits<double>::max();
    extents[1] = extents[3] = extents[5] =  std::numeric_limits<double>::max();

    extents = db->get("init_fission_src", extents);

    // get the low and upper bounds of the geometry
    const cuda_utils::Space_Vector &low_edge  = d_geometry.get_host_ptr()->lower();
    const cuda_utils::Space_Vector &high_edge = d_geometry.get_host_ptr()->upper();

    double lower_x = extents[2 * I];
    double upper_x = extents[2 * I + 1];
    lower_x = std::max(lower_x, low_edge.x);
    upper_x = std::min(upper_x, high_edge.x);
    d_lower.x = lower_x;
    d_width.x = upper_x - lower_x;

    double lower_y = extents[2 * J];
    double upper_y = extents[2 * J + 1];
    lower_y = std::max(lower_y, low_edge.y);
    upper_y = std::min(upper_y, high_edge.y);
    d_lower.y = lower_y;
    d_width.y = upper_y - lower_y;

    double lower_z = extents[2 * K];
    double upper_z = extents[2 * K + 1];
    lower_z = std::max(lower_z, low_edge.z);
    upper_z = std::min(upper_z, high_edge.z);
    d_lower.z = lower_z;
    d_width.z = upper_z - lower_z;

    // store the total number of requested particles per active_cycle
    d_np_requested = static_cast<int>(db->get("Np", 1000));

    // initialize the total for the first active_cycle
    d_np_total = d_np_requested;

    // Allocate device data.
    d_vector_size = db->get("particle_vector_size",10000);
    cuda_utils::memory::Malloc( d_fission_sites_device, d_vector_size );
    cuda_utils::memory::Malloc( d_fission_cells_device, d_vector_size );

    // Get the number of active cycles.
    d_num_inactive_cycle = db->get<int>("num_inactive_cycles");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template <class Geometry>
Fission_Source<Geometry>::~Fission_Source()
{
    cuda_utils::memory::Free( d_fission_sites_device );
    cuda_utils::memory::Free( d_fission_cells_device );
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the initial fission source.
 */
template <class Geometry>
void Fission_Source<Geometry>::build_initial_source()
{
    // send an empty mesh and view
    SDP_Cart_Mesh     mesh;
    Const_Array_View fis_dens;
    build_initial_source(mesh, fis_dens);

    ENSURE(d_wt >= 1.0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the initial source from a mesh distribution.
 */
template <class Geometry>
void Fission_Source<Geometry>::build_initial_source(const SDP_Cart_Mesh& mesh,
                                                    Const_Array_View& fis_dens)
{
    REQUIRE(d_np_total > 0);

    // set the fission site container to an unassigned state
    d_fission_sites = SP_Fission_Sites();
    CHECK(!d_fission_sites);

    // build the domain-replicated fission source
    build_DR(mesh, fis_dens);

    // set counters
    d_np_left = d_np_domain;
    d_np_run  = 0;

    // weight per particle
    d_wt = static_cast<double>(d_np_requested) /
           static_cast<double>(d_np_total);

    // initialize starting cell
    d_current_cell = 0;

    ENSURE(d_wt > 0.0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a source from a fission site container.
 *
 * \param fission_sites
 */
template <class Geometry>
void Fission_Source<Geometry>::build_source(SP_Fission_Sites &fission_sites)
{
    REQUIRE(fission_sites);

    // build an empty fission site container if one doesn't exist (meaning
    // that we haven't yet initialized from an existing fission source)
    if (!d_fission_sites)
    {
        d_fission_sites = std::make_shared<Fission_Site_Container>();
    }

    // the internal fission source should be empty
    REQUIRE(d_fission_sites->empty());

    // swap the input fission sites with the internal storage fission sites
    d_fission_sites.swap(fission_sites);

    // rebalance across sets (when number of blocks per set > 1; the
    // set-rebalance may try to do some load-balancing when it can, that is
    // why this call should comm after the gather; otherwise the
    // load-balancing could provide poor results)
    d_fission_rebalance->rebalance(*d_fission_sites);

    // get the number of fission sites on this domain, on this set, and
    // globally from the fission-rebalance
    d_np_domain = d_fission_rebalance->num_fissions();
    d_np_total  = d_fission_rebalance->num_global_fissions();
    CHECK(d_np_domain >= d_fission_sites->size()); // there could be multiple
                                                    // fissions at a single
                                                    // site

    // set counters
    d_np_left = d_np_domain;
    d_np_run  = 0;
    // weight per particle
    d_wt = static_cast<double>(d_np_requested) /
           static_cast<double>(d_np_total);

    // Increment the active_cycle count.
    ++d_active_cycle_counter;

    ENSURE(d_wt > 0.0);
    ENSURE(fission_sites);
    ENSURE(fission_sites->empty());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a fission site container.
 */
template <class Geometry>
auto Fission_Source<Geometry>::create_fission_site_container() const
    -> SP_Fission_Sites
{
    SP_Fission_Sites fs(std::make_shared<Fission_Site_Container>());
    ENSURE(fs);
    ENSURE(fs->empty());
    return fs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get particles from the source.
*/
template <class Geometry>
void Fission_Source<Geometry>::get_particles(
    cuda_utils::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles )
{
    REQUIRE(d_wt > 0.0);
    REQUIRE( particles.get_host_ptr()->size() == d_vector_size );

    // do nothing if no source
    if (!d_np_left)
    {
        return;
    }

    // Get the particles that are dead.
    int num_particle =
        particles.get_host_ptr()->get_event_size( events::DEAD );

    // Calculate the total number of particles we will create.
    int num_to_create = std::min( d_np_left, num_particle );

    // Get CUDA launch parameters.
    REQUIRE( cuda_utils::Hardware<cuda_utils::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda_utils::Hardware<cuda_utils::arch::Device>::default_block_size();
    unsigned int num_blocks = num_to_create / threads_per_block;
    if ( num_to_create % threads_per_block > 0 ) ++num_blocks;

    // Sample the fission sites if we have them.
    if ( d_fission_sites )
    {
	sample_fission_sites(
	    particles, num_to_create, num_blocks, threads_per_block );
    }

    // If this is an initial source and we have a fission mesh, sample that.
    else if ( d_fis_mesh )
    {
	sample_mesh(
	    particles, num_to_create, num_blocks, threads_per_block );
    }

    // Otherwise this is an initial source and we have no mesh so sample the
    // geometry.
    else
    {
	sample_geometry(
	    particles, num_to_create, num_blocks, threads_per_block );
    }

    // update counters
    d_np_left -= num_to_create;
    d_np_run += num_to_create;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the initial source for Domain Replicated decompositions.
 *
 * In DR decompositions the number of sets is equal to the number of domains
 * (1 block per set).  Thus, the number of particles per set is equal to the
 * number of particles per domain.
 */
template <class Geometry>
void Fission_Source<Geometry>::build_DR(const SDP_Cart_Mesh& mesh,
                                        Const_Array_View& fis_dens)
{
    REQUIRE(profugus::Global_RNG::d_rng.assigned());

    // calculate the number of particles per domain and set (equivalent)
    d_np_domain = d_np_total / profugus::nodes();

    // recalculate the total number of particles (we want the same number of
    // particles in each domain, so the total may change slightly from the
    // requested value)
    d_np_total = d_np_domain * profugus::nodes();

    // if there is a mesh then do stratified sampling to calculate the initial
    // fission distribution
    if (mesh)
    {
        REQUIRE(mesh.get_host_ptr()->num_cells() == fis_dens.size());

        // number of cells in the mesh
        int num_cells = mesh.get_host_ptr()->num_cells();
        // determine the total number of fissions
        double fissions = 0.0;
        for (int cell = 0; cell < num_cells; ++cell)
        {
            fissions += fis_dens[cell] * mesh.get_host_ptr()->volume_host(cell);
        }
        CHECK(fissions > 0.0);
        // allocate fission distribution
        Vec_Int n(num_cells, 0);

        // pre-sample sites on this domain
        double nu            = 0.0;
        int    new_np_domain = 0;
        for (int cell = 0; cell < num_cells; ++cell)
        {
            // calculate the expected number of sites in this cell
            nu = fis_dens[cell] * d_np_domain * 
		 mesh.get_host_ptr()->volume_host(cell) / fissions;

            // there can be n or n+1 sites; with probability n+1-nu there will
            // be n sites, with probability nu-n there will be n+1 sites
            n[cell] = nu;
            if (profugus::Global_RNG::d_rng.ran() < nu - static_cast<double>(n[cell]))
            {
                ++n[cell];
            }

            // add up the number of particles on this domain
            new_np_domain += n[cell];
        }

        // store the distributions persistently
        std::swap(n, d_fis_dist);
        CHECK(d_fis_dist.size() == num_cells);

        // update the number of particles globally and on the domain
        d_np_domain = new_np_domain;
        d_np_total  = d_np_domain;
        profugus::global_sum(&d_np_total, 1);

        // store the mesh for later use
        d_fis_mesh = mesh;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get particles from the fission sites.
 */
template <class Geometry>
void Fission_Source<Geometry>::sample_fission_sites(
    cuda_utils::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles,
    const int num_particle,
    const unsigned int num_blocks,
    const unsigned int threads_per_block )
{
    // Extract the fission sites.
    int copy_start = d_fission_sites->size() - num_particle;
    cuda_utils::memory::Copy_To_Device_Async( d_fission_sites_device,
                                        d_fission_sites->data() + copy_start,
                                        num_particle,
                                        d_stream );

    // Remove the fission sites from the host that we just copied to the
    // device.
    d_fission_sites->resize( copy_start );

    // Launch the kernel.
    sample_fission_sites_kernel<<<num_blocks,threads_per_block,0,d_stream.handle()>>>(
	d_geometry.get_device_ptr(),
	d_physics.get_device_ptr(),
	num_particle,
	d_fission_sites_device,
	d_wt,
        d_active_cycle_counter - d_num_inactive_cycle,
	particles.get_device_ptr() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get particles from the fission mesh.
 */
template <class Geometry>
void Fission_Source<Geometry>::sample_mesh(
    cuda_utils::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles,
    const int num_particle,
    const unsigned int num_blocks,
    const unsigned int threads_per_block )
{
    // Extract and unroll the fission distribution.
    std::vector<int> fission_cells( num_particle );
    int num_extracted = 0;
    int current_cell = 0;
    while ( num_extracted < num_particle )
    {
	CHECK( current_cell < d_fis_dist.size() );

	if ( d_fis_dist[current_cell] )
	{
	    fission_cells[num_extracted] = current_cell;
	    ++num_extracted;
	    --d_fis_dist[current_cell];
	}
	else
	{
	    ++current_cell;
	}
    }

    // Copy the fission cells to the device.
    cuda_utils::memory::Copy_To_Device_Async( d_fission_cells_device,
                                        fission_cells.data(),
                                        num_particle,
                                        d_stream );

    // Launch the kernel.
    sample_mesh_kernel<<<num_blocks,threads_per_block,0,d_stream.handle()>>>(
	d_geometry.get_device_ptr(),
	d_physics.get_device_ptr(),
	d_fis_mesh.get_device_ptr(),
	d_fission_cells_device,
	num_particle,
	d_wt,
        d_active_cycle_counter - d_num_inactive_cycle,
	particles.get_device_ptr() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get particles from the geometry.
 */
template <class Geometry>
void Fission_Source<Geometry>::sample_geometry(
    cuda_utils::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles,
    const int num_particle,
    const unsigned int num_blocks,
    const unsigned int threads_per_block )
{
    // Launch the kernel.
    sample_geometry_kernel<<<num_blocks,threads_per_block,0,d_stream.handle()>>>(
    	d_geometry.get_device_ptr(),
    	d_physics.get_device_ptr(),
    	num_particle,
    	d_wt,
    	d_width,
    	d_lower,
        d_active_cycle_counter - d_num_inactive_cycle,
    	particles.get_device_ptr() );
}

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Fission_Source_t_cuh

//---------------------------------------------------------------------------//
//                 end of Fission_Source.t.cuh
//---------------------------------------------------------------------------//
