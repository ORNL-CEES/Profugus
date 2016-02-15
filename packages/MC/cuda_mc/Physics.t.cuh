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

#include "harness/Soft_Equivalence.hh"
#include "harness/Warnings.hh"

#include "utils/Constants.hh"
#include "utils/Vector_Functions.hh"

#include "cuda_utils/Memory.cuh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Hardware.hh"

#include "Physics.hh"

namespace cuda_profugus
{
//---------------------------------------------------------------------------//
// CUDA KERNELS
//---------------------------------------------------------------------------//
template<class Geometry>
__global__ void initialize_kernel( const double* energy,
				   Particle_Vector<Geometry>* particles )
{
    idx = threadIdx.x + blockIdx.x * blockDim.x;

    // check to make sure the energy is in the group structure and get the
    // group index
    int  group_index = 0;
    bool success     = d_gb.find(energy[idx], group_index);

    // set the group index in the particle
    p.set_group(idx,group_index);
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor that implicitly creates Group_Bounds
 */
template <class Geometry>
Physics<Geometry>::Physics(RCP_Std_DB db,
                           RCP_XS     mat)
    : d_mat(*mat)
    , d_Ng(mat->num_groups())
    , d_Nm(mat->num_mat())
    , d_scatter(d_Nm)
    , d_fissionable(d_Nm)
{
    REQUIRE(!db.is_null());
    REQUIRE(d_mat);
    REQUIRE(d_mat->num_groups() > 0);
    REQUIRE(d_mat->num_mat() > 0);

    // Make the group bounds.
    d_gb = shared_device_ptr<Group_Bounds>(
	Vec_Dbl(mat->bounds().values(),
		mat->bounds().values() + mat->bounds().length()) );
    INSIST(d_gb.num_groups() == mat->num_groups(),
            "Number of groups in material is inconsistent with Group_Bounds.");

    // implicit capture flag
    d_implicit_capture = db->get("implicit_capture", true);

    // check for balanced scattering tables
    d_check_balance = db->get("check_balance", false);

    // turn check balance on if we are not doing implicit capture
    if (!d_implicit_capture) d_check_balance = true;

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
        const auto &sig_s = mat->matrix(matids[m], 0);
        CHECK(sig_s.numRows() == d_Ng);
        CHECK(sig_s.numCols() == d_Ng);

        // loop over all groups and calculate the in-scatter from other
        // groups and add them to the group OUT-SCATTER; remember, we
        // store data as inscatter for the deterministic code
        for (int g = 0; g < d_Ng; g++)
        {
            // get the g column (g->g' scatter stored as g'g in the matrix)
            const auto *column = sig_s[g];

            // add up the scattering
            for (int gp = 0; gp < d_Ng; ++gp)
            {
                matid_scatter[g] += column[gp];
            }
        }

        // check scattering correctness if needed
        if (d_check_balance)
        {n
            for (int g = 0; g < d_Ng; g++)
            {
                if (matid_scatter[g] > mat->vector(matids[m], XS_t::TOTAL)[g])
                {
                    std::ostringstream mm;
                    mm << "Scattering greater than total "
                       << "for material" << m << " in group " << g
                       << ". Total xs is "
                       << mat->vector(matids[m], XS_t::TOTAL)[g]
                       << " and scatter is " << matid_scatter[g];

                    // terminate if we are running analog
                    if (!d_implicit_capture)
                        VALIDATE(false, mm.str());
                    // else add to warnings
                    else
                        ADD_WARNING(mm.str());
                }
            }
        }

	// Copy the scattering to the device.
	offset = m * d_Ng;
	cuda::memory::Copy_To_Device( 
	    d_scatter + offset, matid_scatter.data(), d_Ng );

        // see if this material is fissionable by checking Chi
        host_fissionable[m] = mat->vector(matid, XS_t::CHI).normOne() > 0.0 ?
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
    Shared_Device_Ptr<Particle_Vector_t>& particles )
{
    REQUIRE( energy.size() == particles.get_host_ptr()->size() );

    // Copy the energies to the device.
    double* device_energy;
    cuda::memory::Malloc( device_energy, energy.size() );
    cuda::memory::Copy_To_Device( device_energy, energy.data(), energy.size() );

    // Get CUDA launch parameters.
    REQUIRE( cuda::Hardware<cuda::arch::Device>::have_acquired() );
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = particles.get_host_ptr()->size() / threads_per_block;
    if ( num_particle % threads_per_block > 0 ) ++num_blocks;

    // Initialize the particles.
    initialize_kernel<<<num_blocks, threads_per_block>>>(
	device_energy, particles.get_device_ptr() );

    // Free the device energy array.
    cuda::memory::Free( device_energy );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a total cross section from the physics library.
 */
template <class Geometry>
double Physics<Geometry>::total(physics::Reaction_Type  type,
                                const Particle_t       &p)
{
    REQUIRE(d_mat->num_mat() == d_Nm);
    REQUIRE(d_mat->num_groups() == d_Ng);
    REQUIRE(p.group() < d_Ng);

    // get the matid from the particle
    unsigned int matid = p.matid();
    CHECK(d_mat->has(matid));

    // return the approprate reaction type
    switch (type)
    {
        case physics::TOTAL:
            return d_mat->vector(matid, XS_t::TOTAL)[p.group()];

        case physics::SCATTERING:
            return d_scatter[d_mid2l[matid]][p.group()];

        case physics::FISSION:
            return d_mat->vector(matid, XS_t::SIG_F)[p.group()];

        case physics::NU_FISSION:
            return d_mat->vector(matid, XS_t::NU_SIG_F)[p.group()];

        default:
            return 0.0;
    }

    // undefined or unassigned type, return 0
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a particle through a physical collision.
 */
template <class Geometry>
void Physics<Geometry>::collide(Particle_t &particle)
{
    REQUIRE(d_geometry);
    REQUIRE(particle.event() == events::COLLISION);
    REQUIRE(!d_mat.is_null());
    REQUIRE(d_mat->num_mat() == d_Nm);
    REQUIRE(d_mat->num_groups() == d_Ng);
    REQUIRE(particle.group() < d_Ng);

    // get the material id of the current region
    d_matid = particle.matid();
    CHECK(d_mid2l[static_cast<unsigned int>(d_matid)] < d_Nm);
    CHECK(d_geometry->matid(particle.geo_state()) == d_matid);

    // get the group index
    int group = particle.group();

    // calculate the scattering cross section ratio
    double c = d_scatter[d_mid2l[d_matid]][group] /
               d_mat->vector(d_matid, XS_t::TOTAL)[group];
    CHECK(!d_implicit_capture ? c <= 1.0 : c >= 0.0);

    // we need to do analog transport if the particle is c = 0.0 regardless of
    // whether implicit capture is on or not

    // do implicit capture
    if (d_implicit_capture && c > 0.0)
    {
        // set the event
        particle.set_event(events::IMPLICIT_CAPTURE);

        // do implicit absorption
        particle.multiply_wt(c);
    }

    // do analog transport
    else
    {
        // sample the interaction type
        if (particle.rng().ran() > c)
        {
            // set event indicator
            particle.set_event(events::ABSORPTION);

            // kill particle
            particle.kill();
        }
        else
        {
            // set event indicator
            particle.set_event(events::SCATTER);
        }
    }

    // process scattering events
    if (particle.event() != events::ABSORPTION)
    {
        // determine new group of particle
        group = sample_group(d_matid, group, particle.rng().ran());
        CHECK(group >= 0 && group < d_Ng);

        // set the group
        particle.set_group(group);

        // sample isotropic scattering event
        double costheta = 1.0 - 2.0 * particle.rng().ran();
        double phi      = 2.0 * constants::pi * particle.rng().ran();

        // update the direction of the particle in the geometry-tracker state
        d_geometry->change_direction(costheta, phi, particle.geo_state());
    }
}

//---------------------------------------------------------------------------//
//! Get the energy from a particle via its physics state
template<class Geometry>
double Physics<Geometry>::energy(const Particle_t &p) const
{
    double low = 0.0, up = 0.0;
    d_gb.get_energy(p.group(), low, up);
    return low;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample the fission spectrum and initialize the particle.
 *
 * This function is optimized based on the assumption that nearly all of the
 * fission emission is in the first couple of groups.
 *
 * \post the particle is initialized if fission is sampled
 *
 * \return true if fissionable material and spectrum sampled; false if no
 * fissionable material present
 */
template <class Geometry>
bool Physics<Geometry>::initialize_fission(unsigned int  matid,
                                           Particle_t   &p)
{
    REQUIRE(d_mat->has(matid));

    // sampled flag
    bool sampled = false;

    // only do sampling if this is a fissionable material
    if (is_fissionable(matid))
    {
        // sample the fission group
        int group = sample_fission_group(matid, p.rng().ran());
        sampled   = true;

        // set the group
        p.set_group(group);
    }

    // return the sampled flag
    return sampled;
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
int Physics<Geometry>::sample_fission_site(const Particle_t       &p,
                                           Fission_Site_Container &fsc,
                                           double                  keff)
{
    REQUIRE(d_geometry);
    REQUIRE(d_mat->has(p.matid()));

    // material id
    unsigned int matid = p.matid();

    // if the material is not fissionable exit
    if (!is_fissionable(matid))
        return 0;

    // otherwise make a fission site and sample

    // get the group from the particle
    int group = p.group();

    // calculate the number of fission sites (random number samples to nearest
    // integer)
    int n = static_cast<int>(
        p.wt() *
        d_mat->vector(matid, XS_t::NU_SIG_F)[group] /
        d_mat->vector(matid, XS_t::TOTAL)[group] /
        keff + p.rng().ran());

    // add sites to the fission site container
    for (int i = 0; i < n; ++i)
    {
        Fission_Site site;
        site.m = matid;
        site.r = d_geometry->position(p.geo_state());
        fsc.push_back(site);
    }

    return n;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize a physics state for transport at a previously sampled
 * fission site.
 *
 * This function initializes the physics state at a fission site (so fission
 * will definitely be sampled).  It returns true if the physics state was
 * initialized; it returns false if there are not more particles left to
 * sample at the site.
 *
 * \return true if physics state initialized; false if no particles are left
 * at the site
 */
template <class Geometry>
bool Physics<Geometry>::initialize_fission(Fission_Site &fs,
                                           Particle_t   &p)
{
    REQUIRE(d_mat->has(fs.m));
    REQUIRE(is_fissionable(fs.m));

    // sample the fission group
    int group = sample_fission_group(fs.m, p.rng().ran());
    CHECK(group < d_Ng);

    // set the particle group
    p.set_group(group);

    // we successfully initialized the state
    return true;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Sample a group after a scattering event.
 */
template <class Geometry>
int Physics<Geometry>::sample_group(int    matid,
                                    int    g,
                                    double rnd) const
{
    REQUIRE(!d_mat.is_null());
    REQUIRE(d_mat->num_groups() == d_Ng);
    REQUIRE(g >= 0 && g < d_Ng);
    REQUIRE(rnd >= 0.0 && rnd < 1.0);
    REQUIRE(d_mat->has(matid));
    REQUIRE(d_mid2l.exists(matid));

    // running cdf
    double cdf = 0.0;

    // total out-scattering for this cell and group
    double total = 1.0 / d_scatter[d_mid2l[matid]][g];

    // get the P0 scattering cross section matrix the g column (which is the
    // outscatter) for this group (g->g' is the {A_(g'g) g'=0,Ng} entries of
    // the inscatter matrix
    const auto *scat_g = d_mat->matrix(matid, 0)[g];

    // sample g'
    for (int gp = 0; gp < d_Ng; ++gp)
    {
        // calculate the cdf for scattering to this group
        cdf += scat_g[gp] * total;

        // see if we have sampled this group
        if (rnd <= cdf)
            return gp;
    }
    CHECK(soft_equiv(cdf, 1.0));

    // we failed to sample
    VALIDATE(false, "Failed to sample group.");
    return -1;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a fission group.
 *
 * This function is optimized based on the assumption that nearly all of the
 * fission emission is in the first couple of groups.
 */
template <class Geometry>
int Physics<Geometry>::sample_fission_group(unsigned int matid,
                                            double       rnd) const
{
    REQUIRE(d_mat->has(matid));
    REQUIRE(is_fissionable(matid));

    // running cdf; we make the cdf on the fly because nearly all of the
    // emission is in the first couple of groups so its not worth storing for
    // a binary search
    double cdf = 0.0;

    // get the fission chi
    const auto &chi = d_mat->vector(matid, XS_t::CHI);
    CHECK(chi.length() == d_Ng);

    // sample cdf
    for (int g = 0; g < d_Ng; ++g)
    {
        // update cdf
        cdf += chi[g];

        // check for sampling; update particle's physics state and return
        if (rnd <= cdf)
        {
            // update the group in the particle
            return g;
        }
    }

    // we failed to sample
    VALIDATE(false, "Failed to sample fission group.");
    return -1;
}

} // end namespace cuda_profugus

#endif // cuda_mc_Physics_t_cuh

//---------------------------------------------------------------------------//
//                 end of Physics.t.cuh
//---------------------------------------------------------------------------//
