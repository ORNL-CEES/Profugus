//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Physics.t.cuh
 * \author Stuart Slattery
 * \brief  MG_Physics template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Physics_t_cuh
#define cuda_mc_Physics_t_cuh

#include <sstream>
#include <algorithm>

#include "harness/Warnings.hh"

#include "utils/Constants.hh"
#include "utils/Vector_Functions.hh"

#include "cuda_utils/Memory.cuh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Hardware.hh"
#include "cuda_utils/Utility_Functions.hh"

#include "Physics.hh"

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
// CUDA DEVICE FUNCTIONS
//---------------------------------------------------------------------------//
/*
 * \brief Sample a group.
 */
__device__ int sample_group( const XS_Device* xs,
			     const double* scatter,
			     const int* matid_g2l,
			     const int matid,
			     const int g,
			     const double rnd )
{
    // running cdf
    double cdf = 0.0;

    // total out-scattering for this cell and group
    int num_groups = xs->num_groups();
    double total = 1.0 / scatter[matid_g2l[matid] * num_groups + g];

    // get the P0 scattering cross section matrix the g column (which is the
    // outscatter) for this group (g->g' is the {A_(g'g) g'=0,Ng} entries of
    // the inscatter matrix
    const auto scat_matrix = xs->matrix(matid, 0);

    // sample g'
    for (int gp = 0; gp < num_groups; ++gp)
    {
        // calculate the cdf for scattering to this group
        cdf += scat_matrix(gp,g) * total;

        // see if we have sampled this group
        if (rnd <= cdf)
            return gp;
    }
    CHECK( cuda::utility::soft_equiv(cdf, 1.0) );

    // we failed to sample
    return -1;
}

//---------------------------------------------------------------------------//
// CUDA GLOBAL KERNELS
//---------------------------------------------------------------------------//
/*
 * \brief initialize particles with a given energy
 */
template<class Geometry>
__global__ void initialize_kernel( const double* energy,
				   const Group_Bounds* gb,
				   Particle_Vector<Geometry>* particles )
{
    // get the thread id
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    // this kernel updates all particles
    if ( idx < particles->size() )
    {
	// check to make sure the energy is in the group structure and get the
	// group index
	int  group_index = 0;
	bool success = gb->find( energy[idx], group_index );

	// set the group index in the particle
	particles->set_group( idx, group_index );
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Process particles through a collision.
 */
template<class Geometry>
__global__ void collide_kernel( const std::size_t start_idx,
				const std::size_t num_particle,
				const Geometry* geometry,
				const XS_Device* xs,
				const int* matid_g2l,
				const double* scatter,
				const bool implicit_capture,
				Particle_Vector<Geometry>* particles )
{
    // get the thread id.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if ( idx < num_particle )
    {
	// get the particle index
	int pidx = idx + start_idx;
	
	REQUIRE(geometry);
	REQUIRE(particles->event(pidx) == events::COLLISION);

	// get the material id of the current region
	int matid = particles->matid(pidx);
	CHECK(geometry->matid(particles->geo_state(pidx)) == matid);

	// get the group index
	int group = particles->group(pidx);

	// calculate the scattering cross section ratio
	double c = scatter[matid_g2l[matid]*xs->num_groups() + group] /
		   xs->vector(matid, profugus::XS::TOTAL)(group);
	CHECK(!implicit_capture ? c <= 1.0 : c >= 0.0);

	// we need to do analog transport if the particles->is c = 0.0
	// regardless of whether implicit capture is on or not

	// do implicit capture
	if (implicit_capture && c > 0.0)
	{
	    // set the event
	    particles->set_event(pidx,events::IMPLICIT_CAPTURE);

	    // do implicit absorption
	    particles->multiply_wt(pidx,c);
	}

	// do analog transport
	else
	{
	    // sample the interaction type
	    if (particles->ran(pidx) > c)
	    {
		// set event indicator
		particles->set_event(pidx,events::ABSORPTION);

		// kill particle
		particles->kill(pidx);
	    }
	    else
	    {
		// set event indicator
		particles->set_event(pidx,events::SCATTER);
	    }
	}

	// process scattering events
	if (particles->event(pidx) != events::ABSORPTION)
	{
	    // determine new group of particle
	    group = sample_group(xs, scatter, matid_g2l,
				 matid, group, particles->ran(pidx));
	    CHECK(group >= 0 && group < xs->num_groups());

	    // set the group
	    particles->set_group(pidx,group);

	    // sample isotropic scattering event
	    double costheta = 1.0 - 2.0 * particles->ran(pidx);
	    double phi      = 4.0 * std::asin(1.0) * particles->ran(pidx);

	    // update the direction of the particles->in the geometry-tracker
	    // state
	    geometry->change_direction(
		costheta, phi, particles->geo_state(pidx));
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample fission.
 */
template<class Geometry>
__global__ void sample_fission_site_kernel(
    const std::size_t start_idx,
    const std::size_t num_particle,
    const Geometry* geometry,
    const XS_Device* xs,
    const double keff,
    const int* matid_g2l,    
    const int* is_fissionable,
    Particle_Vector<Geometry>* particles,
    typename Physics<Geometry>::Device_Fission_Site* fission_sites )
{
    // get the thread id.
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if ( idx < num_particle )
    {
	// get the particle index
	int pidx = idx + start_idx;

	// material id
	unsigned int matid = particles->matid(pidx);

	typename Physics<Geometry>::Device_Fission_Site& site =
	    fission_sites[idx];

	// if fissionable calculate the number of fission sites
	if ( is_fissionable[matid_g2l[matid]] )
	{
	    // get the group from the particle
	    int group = particles->group(pidx);

	    // calculate the number of fission sites (random number samples to
	    // nearest integer)
	    int n = static_cast<int>(
		particles->wt(pidx) *
		xs->vector(matid, profugus::XS::NU_SIG_F)(group) /
		xs->vector(matid, profugus::XS::TOTAL)(group) /
		keff + particles->ran(pidx) );

	    // add sites to the fission site container
	    site.m = matid;
	    site.r = geometry->position(particles->geo_state(pidx));
	    site.n = n;
	}

	// otherwise there is no fission
	else
	{
	    site.n = 0;
	}
    }
}					    

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor that implicitly creates Group_Bounds
 */
template <class Geometry>
Physics<Geometry>::Physics( ParameterList_t& db, const XS_t& mat )
    : d_mat_device( nullptr )
    , d_Ng(mat.num_groups())
    , d_Nm(mat.num_mat())
    , d_gb_device( nullptr )
    , d_matid_g2l( nullptr )
    , d_scatter( nullptr )
    , d_fissionable( nullptr )
    , d_device_sites( nullptr )     
{
    // Create the device cross sections
    d_mat = cuda::shared_device_ptr<XS_Device>( mat );
    d_mat_device = d_mat.get_device_ptr();
    
    REQUIRE(d_mat);
    REQUIRE(mat.num_groups() > 0);
    REQUIRE(mat.num_mat() > 0);

    // Make the group bounds.
    d_gb = cuda::shared_device_ptr<Group_Bounds>(
	Vec_Dbl(mat.bounds().values(),
		mat.bounds().values() + mat.bounds().length()) );
    d_gb_device = d_gb.get_device_ptr();
    INSIST(d_gb.get_host_ptr()->num_groups() == mat.num_groups(),
            "Number of groups in material is inconsistent with Group_Bounds.");

    // implicit capture flag
    d_implicit_capture = db.get("implicit_capture", true);

    // Get the matids.
    Vec_Int matids;
    mat.get_matids( matids );

    // Create a global to local mapping of matids.
    int matid_g2l_size = *std::max_element( matids.begin(), matids.end() ) + 1;
    Teuchos::Array<int> host_matid_g2l( matid_g2l_size, -1 );
    for ( int m = 0; m < d_Nm; ++m )
    {
	host_matid_g2l[ matids[m] ] = m;
    }

    // Allocate a matid global-to-local map.
    cuda::memory::Malloc( d_matid_g2l, matid_g2l_size );

    // Copy the matid list to the device.
    cuda::memory::Copy_To_Device( 
	d_matid_g2l, host_matid_g2l.getRawPtr(), matid_g2l_size );
    host_matid_g2l.clear();

    // Allocate scattering.
    cuda::memory::Malloc( d_scatter, d_Nm * d_Ng );

    // calculate total scattering over all groups for each material and
    // determine if fission is available for a given material
    std::size_t offset = 0;
    std::vector<double> matid_scatter( d_Ng, 0.0 );
    std::vector<int> host_fissionable( d_Nm, false );
    for ( int m = 0; m < d_Nm; ++m )
    {
	// Clear the scatter vector.
	matid_scatter.assign( d_Ng, 0.0 );

        // get the P0 scattering matrix for this material
        const auto &sig_s = mat.matrix(matids[m], 0);
        CHECK(sig_s.numRows() == d_Ng);
        CHECK(sig_s.numCols() == d_Ng);

        // loop over all groups and calculate the in-scatter from other
        // groups and add them to the group OUT-SCATTER; remember, we
        // store data as inscatter for the deterministic code
        for (int g = 0; g < d_Ng; g++)
        {
            // add up the scattering
            for (int gp = 0; gp < d_Ng; ++gp)
            {
                matid_scatter[g] += sig_s(gp,g);
            }
        }

	// Copy the scattering to the device.
	offset = m * d_Ng;
	cuda::memory::Copy_To_Device( 
	    d_scatter + offset, matid_scatter.data(), d_Ng );

        // see if this material is fissionable by checking Chi
        host_fissionable[m] = mat.vector(m, XS_t::CHI).normOne() > 0.0 ?
			      true : false;
    }

    // Copy the fissionable vector to the device.
    cuda::memory::Malloc( d_fissionable, d_Nm );
    cuda::memory::Copy_To_Device( 
	d_fissionable, host_fissionable.data(), d_Nm );

    ENSURE(d_Nm > 0);
    ENSURE(d_Ng > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor
 */
template <class Geometry>
Physics<Geometry>::~Physics()
{
    cuda::memory::Free( d_matid_g2l );
    cuda::memory::Free( d_scatter );
    cuda::memory::Free( d_fissionable );
    if ( nullptr != d_device_sites )
    {
	cuda::memory::Free( d_device_sites );
    }
}

//---------------------------------------------------------------------------//
// DERIVED PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize the physics state.
 *
 * The correct group corresponding to E is determined for a particle.
 *
 * \param energy energy in eV
 * \param p particle
 */
template <class Geometry>
void Physics<Geometry>::initialize(
    const std::vector<double>& energy, 
    cuda::Shared_Device_Ptr<Particle_Vector_t>& particles )
{
    int num_particle = particles.get_host_ptr()->size();
    REQUIRE( energy.size() == num_particle );

    // Copy the energies to the device.
    double* device_energy;
    cuda::memory::Malloc( device_energy, energy.size() );
    cuda::memory::Copy_To_Device( device_energy, energy.data(), energy.size() );

    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = num_particle / threads_per_block;
    if ( num_particle % threads_per_block > 0 ) ++num_blocks;

    // Initialize the particles.
    initialize_kernel<<<num_blocks, threads_per_block>>>(
	device_energy, d_gb.get_device_ptr(), particles.get_device_ptr() );

    // Free the device energy array.
    cuda::memory::Free( device_energy );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a particle through a physical collision.
 */
template <class Geometry>
void Physics<Geometry>::collide(
    cuda::Shared_Device_Ptr<Particle_Vector_t>& particles )
{
    // Get the particles that will have a collision.
    std::size_t start_idx = 0;
    std::size_t num_particle = 0;
    particles.get_host_ptr()->get_event_particles( events::COLLISION,
						   start_idx,
						   num_particle );
    
    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = num_particle / threads_per_block;
    if ( num_particle % threads_per_block > 0 ) ++num_blocks;

    // Process the collisions.
    collide_kernel<<<num_particle,threads_per_block>>>(
	start_idx,
	num_particle,
	d_geometry.get_device_ptr(),
	d_mat.get_device_ptr(),
	d_matid_g2l,
	d_scatter,
	d_implicit_capture,
	particles.get_device_ptr() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a fission site.
 *
 * A fission site is sampled using the MCNP5 routine:
 * \f[
 n = w\bar{\nu}\frac{\sigma_{\mathrm{f}}}{\sigma_{\mathrm{t}}}
 \frac{1}{k_{\mathrm{eff}}} + \xi\:,
 * \f]
 * where
 * \f[
 \begin{array}{lll}
 w &=& \mbox{weight before analog or implicit capture}\\
 \bar{\nu} &=& \mbox{average neutrons produced by fission}\\
 \sigma_{\mathrm{f}} &=& \mbox{fission cross section}\\
 \sigma_{\mathrm{t}} &=& \mbox{total cross section}\\
 k_{\mathrm{eff}} &=& \mbox{latest eigenvalue iterate}
 \end{array}
 * \f]
 * and \e n is the number of fission events at the site rounded to the nearest
 * integer.
 *
 * \return the number of fission events added at the site
 */
template <class Geometry>
int Physics<Geometry>::sample_fission_site(
    cuda::Shared_Device_Ptr<Particle_Vector_t>& particles,
    Fission_Site_Container &fsc,
    double                  keff)
{
    // Lazy allocate the fission site work vectors.
    if ( nullptr == d_device_sites )
    {
	d_host_sites.resize( particles.get_host_ptr()->size() );
	cuda::memory::Malloc( d_device_sites, particles.get_host_ptr()->size() );
    }

    // Get the particles that will have a collision.
    std::size_t start_idx = 0;
    std::size_t num_particle = 0;
    particles.get_host_ptr()->get_event_particles( events::COLLISION,
						   start_idx,
						   num_particle );

    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = num_particle / threads_per_block;
    if ( num_particle % threads_per_block > 0 ) ++num_blocks;

    // Sample the fission sites.
    sample_fission_site_kernel<<<num_particle,threads_per_block>>>(
	start_idx,
	num_particle,
	d_geometry.get_device_ptr(),
	d_mat.get_device_ptr(),
	keff,
	d_matid_g2l,
	d_fissionable,
	particles.get_device_ptr(),
	d_device_sites );

    // Pull the fission sites off of the device.
    cuda::memory::Copy_To_Host(
	d_host_sites.data(), d_device_sites, num_particle );

    // Add the fission sites to the fission container.
    int num_sites = 0;
    for ( int p = 0; p < num_particle; ++p )
    {
	for ( int n = 0; n < d_host_sites[p].n; ++n )
	{
	    Fission_Site site;
	    site.m = d_host_sites[p].m;
	    site.r = d_host_sites[p].r;
	    fsc.push_back( site );
	}
	
	num_sites += d_host_sites[p].n;
    }

    return num_sites;
}

//---------------------------------------------------------------------------//

} // End namespace cuda_profugus

#endif // cuda_mc_Physics_t_cuh

//---------------------------------------------------------------------------//
//                 end of Physics.t.cuh
//---------------------------------------------------------------------------//
