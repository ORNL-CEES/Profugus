//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Physics.i.cuh
 * \author Thomas M. Evans
 * \date   Thursday May 1 11:14:55 2014
 * \brief  Physics template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Physics_i_cuh
#define cuda_mc_Physics_i_cuh

#include <sstream>
#include <algorithm>

#include "harness/Soft_Equivalence.hh"
#include "harness/Warnings.hh"
#include "utils/Constants.hh"
#include "utils/Vector_Functions.hh"
#include "cuda_utils/Constants.hh"
#include "Physics.cuh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
/*!
 * \brief Process a particle through a physical collision.
 */
template <class Geometry>
__device__
void Physics<Geometry>::collide(int pid, Particle_Vector_t &particles) const
{
    DEVICE_REQUIRE(d_geometry);
    DEVICE_REQUIRE(particles.event(pid) == profugus::events::COLLISION);
    DEVICE_REQUIRE(d_mat);
    DEVICE_REQUIRE(d_mat->num_mat() == d_Nm);
    DEVICE_REQUIRE(d_mat->num_groups() == d_Ng);
    DEVICE_REQUIRE(particles.group(pid) < d_Ng);

    // get the material id of the current region
    int matid = particles.matid(pid);
    DEVICE_CHECK(matid < d_Nm);
    DEVICE_CHECK(d_geometry->matid(particles.geo_state(pid)) == matid);

    // get the group index
    int group = particles.group(pid);

    // calculate the scattering cross section ratio
    double c = d_scatter[group_mat_index(group,matid)] /
               d_mat->vector(matid, XS_t::TOTAL)(group);
    DEVICE_CHECK(!d_implicit_capture ? c <= 1.0 : c >= 0.0);

    // we need to do analog transport if the particle is c = 0.0 regardless of
    // whether implicit capture is on or not

    // do implicit capture
    if (d_implicit_capture && c > 0.0)
    {
        // set the event
        particles.set_event(pid,profugus::events::IMPLICIT_CAPTURE);

        // do implicit absorption
        particles.multiply_wt(pid,c);
    }


    // do analog transport
    else
    {
        // sample the interaction type
        if (particles.ran(pid) > c)
        {
            // set event indicator
            particles.set_event(pid,profugus::events::ABSORPTION);

            // kill particle
            particles.kill(pid);
        }
        else
        {
            // set event indicator
            particles.set_event(pid,profugus::events::SCATTER);
        }
    }

    // process scattering events
    if (particles.event(pid) != profugus::events::ABSORPTION)
    {
        // determine new group of particle
        group = sample_group(matid, group, particles.ran(pid));
        DEVICE_CHECK(group >= 0 && group < d_Ng);

        // set the group
        particles.set_group(pid,group);

        // sample isotropic scattering event
        double costheta = 1.0 - 2.0 * particles.ran(pid);
        double phi      = 2.0 * cuda::constants::pi * particles.ran(pid);

        // update the direction of the particle in the geometry-tracker state
        d_geometry->change_direction(costheta, phi, particles.geo_state(pid));
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a total cross section from the physics library.
 */
template <class Geometry>
__device__
double Physics<Geometry>::total(Reaction_Type            type,
                                int                      pid,
                                const Particle_Vector_t &particles) const
{
    DEVICE_REQUIRE(d_mat->num_mat() == d_Nm);
    DEVICE_REQUIRE(d_mat->num_groups() == d_Ng);
    DEVICE_REQUIRE(particles.group(pid) < d_Ng);

    // get the matid from the particle
    unsigned int matid = particles.matid(pid);
    DEVICE_CHECK( matid < d_mat->num_mat() );

    // return the approprate reaction type
    switch (type)
    {
        case profugus::physics::TOTAL:
            return d_mat->vector(matid, XS_t::TOTAL)(particles.group(pid));

        case profugus::physics::SCATTERING:
            return d_scatter[group_mat_index(particles.group(pid),matid)];

        case profugus::physics::FISSION:
            return d_mat->vector(matid, XS_t::SIG_F)(particles.group(pid));

        case profugus::physics::NU_FISSION:
            return d_mat->vector(matid, XS_t::NU_SIG_F)(particles.group(pid));

        default:
            return 0.0;
    }
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
__device__
bool Physics<Geometry>::initialize_fission(unsigned int       matid,
                                           int                pid,
                                           Particle_Vector_t &particles) const
{
    DEVICE_REQUIRE( matid < d_mat->num_mat() );

    // sampled flag
    bool sampled = false;

    // only do sampling if this is a fissionable material
    if (is_fissionable(matid))
    {
        // sample the fission group
        int group = sample_fission_group(matid,particles.ran(pid));
        sampled   = true;

        // set the group
        particles.set_group(pid,group);
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
__device__
int Physics<Geometry>::sample_fission_site(int                      pid,
                                           const Particle_Vector_t &particles,
                                           double                   keff) const
{
    DEVICE_REQUIRE(d_geometry);
    DEVICE_REQUIRE(particles.matid(pid) < d_mat->num_mat() );

    // material id
    unsigned int matid = particles.matid(pid);

    // if the material is not fissionable exit
    if (!is_fissionable(matid))
        return 0;

    // get the group from the particle
    int group = particles.group(pid);

    // calculate the number of fission sites (random number samples to nearest
    // integer)
    int n = static_cast<int>(
        particles.wt(pid) *
        d_mat->vector(matid, XS_t::NU_SIG_F)(group) /
        d_mat->vector(matid, XS_t::TOTAL)(group) /
        keff + particles.ran(pid));

    // NOTE: Unlike cpu version, we put creation of actual fission site
    // into Domain_Transporter

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
__device__
bool Physics<Geometry>::initialize_fission(const Fission_Site &fs,
                                           int                 pid,
                                           Particle_Vector_t  &particles) const
{
    DEVICE_REQUIRE(fs.m < d_mat->num_mat());
    DEVICE_REQUIRE(is_fissionable(fs.m));

    // sample the fission group
    int group = sample_fission_group(fs.m, particles.ran(pid));
    DEVICE_CHECK(group < d_Ng);

    // set the particle group
    particles.set_group(pid,group);

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
__device__
int Physics<Geometry>::sample_group(int    matid,
                                    int    g,
                                    double rnd) const
{
    DEVICE_REQUIRE(d_mat);
    DEVICE_REQUIRE(d_mat->num_groups() == d_Ng);
    DEVICE_REQUIRE(g >= 0 && g < d_Ng);
    DEVICE_REQUIRE(rnd >= 0.0 && rnd < 1.0);
    DEVICE_REQUIRE(matid < d_mat->num_mat());

    // running cdf
    double cdf = 0.0;

    // total out-scattering for this cell and group
    double total = 1.0 / d_scatter[group_mat_index(g,matid)];

    // get the P0 scattering cross section matrix the g column (which is the
    // outscatter) for this group (g->g' is the {A_(g'g) g'=0,Ng} entries of
    // the inscatter matrix
    const auto & scat = d_mat->matrix(matid,0);

    // sample g'
    for (int gp = 0; gp < d_Ng; ++gp)
    {
        // calculate the cdf for scattering to this group
        cdf += scat(gp,g) * total;

        // see if we have sampled this group
        if (rnd <= cdf)
            return gp;
    }

    // we failed to sample
    DEVICE_INSIST(false, "Failed to sample group.");
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
__device__
int Physics<Geometry>::sample_fission_group(unsigned int matid,
                                            double       rnd) const
{
    DEVICE_REQUIRE(matid < d_mat->num_mat());
    DEVICE_REQUIRE(is_fissionable(matid));

    // running cdf; we make the cdf on the fly because nearly all of the
    // emission is in the first couple of groups so its not worth storing for
    // a binary search
    double cdf = 0.0;

    // get the fission chi
    auto chi = d_mat->vector(matid, XS_t::CHI);
    DEVICE_CHECK(chi.size() == d_Ng);

    // sample cdf
    for (int g = 0; g < d_Ng; ++g)
    {
        // update cdf
        cdf += chi(g);

        // check for sampling; update particle's physics state and return
        if (rnd <= cdf)
        {
            // update the group in the particle
            return g;
        }
    }

    // we failed to sample
    DEVICE_INSIST(false, "Failed to sample fission group.");
    return -1;
}

} // end namespace cuda_mc

#endif // cuda_mc_Physics_i_cuh

//---------------------------------------------------------------------------//
//                 end of Physics.i.cuh
//---------------------------------------------------------------------------//
