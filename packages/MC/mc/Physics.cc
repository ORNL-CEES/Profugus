//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Physics.cc
 * \author Thomas M. Evans
 * \date   Tue Jan 04 12:50:33 2011
 * \brief  MG_Physics member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Physics_cc
#define mc_Physics_cc

#include <sstream>
#include <algorithm>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "harness/Warnings.hh"
#include "utils/Constants.hh"
#include "mc/Definitions.hh"
#include "xslib/AMPX_MT.hh"
#include "utils/Vector_Functions.hh"
#include "geometry/Definitions.hh"
#include "MG_Physics.hh"

namespace shift
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor that implicitly creates Group_Bounds
 */
template<class Geometry>
MG_Physics<Geometry>::MG_Physics(SP_Std_DB db,
                                 SP_XS_DB  mat)
{
    Require (db);
    Require (mat);

    SP_Group_Bounds gb = mc::Group_Bounds::from_db(db);
    Check(gb);

    construct(db, mat, gb);

    Ensure (d_Nm > 0);
    Ensure (d_Ng > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor that implicitly creates Group_Bounds
 */
template<class Geometry>
MG_Physics<Geometry>::MG_Physics(SP_Std_DB       db,
                                 SP_XS_DB        mat,
                                 SP_Group_Bounds gb)
{
    Require (db);
    Require (mat);
    Require (gb);

    construct(db, mat, gb);

    Ensure (d_Nm > 0);
    Ensure (d_Ng > 0);
}

//---------------------------------------------------------------------------//

template<class Geometry>
void MG_Physics<Geometry>::construct(SP_Std_DB       db,
                                     SP_XS_DB        mat,
                                     SP_Group_Bounds gb)
{
    Require (db);
    Require (mat);
    Require (gb);
    Require (mat->num_groups() > 0);
    Require (mat->num_mat() > 0);

    // Assign groups, materials, material numbers
    d_mat = mat;
    d_gb  = gb;
    d_Ng  = d_mat->num_groups();
    d_Nm  = d_mat->num_mat();

    Insist (gb->num_groups() == mat->num_groups(),
            "Number of groups in material is inconsistent with Group_Bounds.");

    // implicit capture flag
    db->default_key("implicit_capture", true);
    d_implicit_capture = db->template get<bool>("implicit_capture");

    // check for balanced scattering tables
    db->default_key("check_balance", false);
    d_check_balance = db->template get<bool>("check_balance");

    // turn check balance on if we are not doing implicit capture
    if (!d_implicit_capture) d_check_balance = true;

    // resize the total group scattering table
    d_scatter.resize(d_Nm);

    // calculate total scattering over all groups for each material
    for (int m = 0; m < d_Nm; ++m)
    {
        // only calculate for materials that are assigned
        if (d_mat->assigned(m))
        {
            // size the group vector for this material
            d_scatter[m].resize(d_Ng, 0.0);

            // loop over all groups and calculate the in-scatter from other
            // groups and add them to the group OUT-SCATTER; remember, we
            // store data as inscatter for the deterministic code
            for (int g = 0; g < d_Ng; g++)
            {
                // downscattering inscatter cross sections for this group
                const typename XS_DB_t::XS &xs = d_mat->get(m, g);

                // process the down-scattering cross sections
                for (int gp = 0; gp < g+1; gp++)
                {
                    // get the inscatter (and within-group) scatter cross
                    // section and add it to the groups outscatter
                    d_scatter[m][gp] += xs.sigma_s(gp, 0);
                }

                // process the upscatter cross sections
                const Vec_Int &columns = xs.upscatter_columns();
                Check (columns.size() == xs.num_upscatter());
                for (int n = 0, N = columns.size(); n < N; n++)
                {
                    // gp -> columns[n] for upscatter we only loop over the
                    // columns that actually have upscatter groups
                    d_scatter[m][columns[n]] += xs.sigma_s(columns[n], 0);
                }
            }

            // check scattering correctness if needed
            if (d_check_balance)
            {
                for (int g = 0; g < d_Ng; g++)
                {
                    if (d_scatter[m][g] > d_mat->get(m, g).sigma())
                    {
                        std::ostringstream mm;
                        mm << "Scattering greater than total "
                           << "for material" << m << " in group " << g
                           << ". Total xs is " << d_mat->get(m,g).sigma()
                           << " and scatter is " << d_scatter[m][g];

                        // terminate if we are running analog
                        if (!d_implicit_capture)
                            Validate(false, mm.str());
                        // else add to warnings
                        else
                            ADD_WARNING(mm.str());
                    }
                }
            }
        }
    }

    // make an transpose graph so that the cross sections are for out-scatter
    // instead of in-scatter
    d_graph.construct(*d_mat);
    d_graph.transpose();

    Ensure (d_graph.adjoint());
}

//---------------------------------------------------------------------------//
// DERIVED PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize the physics state.
 *
 * The correct group corresponding to E is determined for a particle.
 *
 * \param type particle type (eg. mc::NEUTRON)
 * \param energy energy in eV
 * \param state physics state (MG_State)
 */
template<class Geometry>
void MG_Physics<Geometry>::initialize(mc::Particle_Type type,
                                      double            energy,
                                      Physics_State_t&  state)
{
    Require (type < mc::END_PARTICLE_TYPE);

    // set the particle type in the state
    state.type  = type;
    bool success = d_gb->find(type, energy, state.group);

    Validate(success, "Particle " << type << " with energy " << energy
            << " is outside the multigroup boundaries.");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a particle through a physical collision.
 */
template<class Geometry>
void MG_Physics<Geometry>::collide(Particle_t &particle,
                                   Bank_t     &bank)
{
    Require (b_geometry);
    Require (particle.event() == mc::events::COLLISION);
    Require (d_mat);
    Require (d_mat->num_mat() == d_Nm);
    Require (d_mat->num_groups() == d_Ng);

    // return if the particle is not at a collision site

    // get a reference to the particle state
    Physics_State_t &state = particle.physics_state();
    Check (state.group < d_Ng);

    // get the material id of the current region
    d_matid = particle.matid();
    Check (d_matid < d_Nm);
    Check (b_geometry->matid(particle.geo_state()) == d_matid);

    // calculate the scattering cross section ratio
    register double c = d_scatter[d_matid][state.group] /
                        d_mat->get(d_matid, state.group).sigma();
    Check (!d_implicit_capture ? c <= 1.0 : c >= 0.0);

    // we need to do analog transport if the particle is c = 0.0 regardless of
    // whether implicit capture is on or not

    // do implicit capture
    if (d_implicit_capture && c > 0.0)
    {
        // set the event
        particle.set_event(mc::events::IMPLICIT_CAPTURE);

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
            particle.set_event(mc::events::ABSORPTION);

            // kill particle
            particle.kill();
        }
        else
        {
            // set event indicator
            particle.set_event(mc::events::SCATTER);
        }
    }

    // process scattering events
    if (particle.event() != mc::events::ABSORPTION)
    {
        // determine new group of particle
        state.group = sample_group(d_matid, state.group, particle.rng().ran());
        Check (state.group >= 0 && state.group < d_Ng);

        // sample isotropic scattering event
        double costheta = 1.0 - 2.0 * particle.rng().ran();
        double phi      = 2.0 * nemesis::constants::pi * particle.rng().ran();

        // update the direction of the particle in the geometry-tracker state
        b_geometry->change_direction(costheta, phi, particle.geo_state());
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a total cross section from the physics library.
 */
template<class Geometry>
double MG_Physics<Geometry>::total(mc::physics::Reaction_Type  type,
                                   int                         matid,
                                   const Physics_State_t       &state)
{
    Require (d_mat->num_mat() == d_Nm);
    Require (d_mat->num_groups() == d_Ng);
    Require (matid < d_Nm);
    Require (state.group < d_Ng);

    // return the approprate reaction type
    switch (type)
    {
        case mc::physics::TOTAL:
            return d_mat->get(matid, state.group).sigma();

        case mc::physics::SCATTERING:
            return d_scatter[matid][state.group];

        case mc::physics::FISSION:
            // zero if non-fission material
            if (!d_mat->assigned_fission(matid))
                return 0.;

            return d_mat->get_fission(matid).data(
                    state.group,
                    denovo::Fission_Data::SIGMA_F);

        case mc::physics::NU_FISSION:
            // zero if non-fission material
            if (!d_mat->assigned_fission(matid))
                return 0.;

            return d_mat->get_fission(matid).data(
                    state.group,
                    denovo::Fission_Data::NU_SIGMA_F);

        case mc::physics::KAPPA_SIGMA:
            // zero if non-assigned
            if (!d_mat->assigned_xs_opt(matid))
                return 0.;

            return d_mat->get_xs_opt(matid).data(
                    state.group,
                    denovo::Optional_XS_Data::KAPPA_SIGMA);

        default:
            return 0.;
    }

    // undefined or unassigned type, return 0
    return 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample the fission spectrum and initialize the physics state.
 *
 * This function is optimized based on the assumption that nearly all of the
 * fission emission is in the first couple of groups.
 *
 * \post the physics state is initialized if fission is sampled
 *
 * \return true if fissionable material and spectrum sampled; false if no
 * fissionable material present
 */
template<class Geometry>
bool MG_Physics<Geometry>::initialize_fission(int              matid,
                                              RNG              rng,
                                              Physics_State_t &state)
{
    Require (matid < d_mat->num_mat());

    // sampled flag
    bool sampled = false;

    // only do sampling if this is a fissionable material
    if (d_mat->assigned_fission(matid))
    {
        // sample the fission group
        state.group = sample_fission_group(matid, rng.ran());
        state.type  = mc::NEUTRON;
        sampled     = true;
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
template<class Geometry>
int MG_Physics<Geometry>::sample_fission_site(const Particle_t       &particle,
                                              Fission_Site_Container &fsc,
                                              double                  keff)
{
    Require (b_geometry);
    Require (particle.matid() < d_mat->num_mat());

    // material id
    int matid = particle.matid();

    // if the material is not fissionable exit
    if (!d_mat->assigned_fission(matid))
        return 0;

    // otherwise make a fission site and sample

    // get the physics state handle from the particle
    const Physics_State_t &s = particle.physics_state();

    // calculate the number of fission sites (random number samples to nearest
    // integer)
    int n = static_cast<int>(
        particle.wt() *
        d_mat->get_fission(matid).data(s.group, XS_DB_t::FIS_XS::NU_SIGMA_F) /
        d_mat->get(matid, s.group).sigma() /
        keff +
        particle.rng().ran());

    // add sites to the fission site container
    for (int i = 0; i < n; ++i)
    {
        Fission_Site site;
        site.m = matid;
        site.r = b_geometry->position(particle.geo_state());
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
template<class Geometry>
bool MG_Physics<Geometry>::initialize_fission(Fission_Site    &fs,
                                              RNG              rng,
                                              Physics_State_t &state)
{
    Require (fs.m < d_mat->num_mat());
    Require (d_mat->assigned_fission(fs.m));

    // sample the fission group
    state.group = sample_fission_group(fs.m, rng.ran());

    // set the particle type
    state.type = mc::NEUTRON;

    // we successfully initialized the state
    return true;
}
//---------------------------------------------------------------------------//
/*!
 * brief Get minimum energy for given particle type
 *
 * \param type Particle type
 */
template<class Geometry>
double MG_Physics<Geometry>::min_energy(mc::Particle_Type type) const
{
    double energy = 0.0;
    bool success = d_gb->has_particle(type);
    if (success) {
        double num_groups = d_gb->num_groups(type);
        energy = d_gb->group_bounds(type)[num_groups];
    }
    else
        Validate(success, "Particle type not in Group Bounds: " << type);

    return energy;
}

//---------------------------------------------------------------------------//
/*!
 * brief Get maximum energy for given particle type
 *
 * \param type Particle type
 */
template<class Geometry>
double MG_Physics<Geometry>::max_energy(mc::Particle_Type type) const
{
    double energy = 0.0;
    bool success = d_gb->has_particle(type);
    if (success)
        energy = d_gb->group_bounds(type)[0];
    else
        Validate(success, "Particle type not in Group Bounds: " << type);

    return energy;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Sample a group after a scattering event.
 */
template<class Geometry>
int MG_Physics<Geometry>::sample_group(int    matid,
                                       int    g,
                                       double rnd) const
{
    Require (d_mat);
    Require (d_mat->num_groups() == d_Ng);
    Require (g >= 0 && g < d_Ng);
    Require (rnd >= 0.0 && rnd <= 1.0);
    Require (d_mat->assigned(matid));
    Require (matid < d_scatter.size());

    // running cdf
    double cdf = 0.0;

    // total out-scattering for this cell and group
    double total = 1.0 / d_scatter[matid][g];

    // don't allow hard 1.0
    if (rnd == 1.0) rnd -= 1.0e-12;

    // determine the upscatter out-scattering columns for this group (from the
    // graph)
    const Vec_Int &columns = d_graph.lower(g);

    // check upscatter columns first because they start at column 0
    for (int n = 0, N = columns.size(); n < N; ++n)
    {
        // cross sections
        const typename XS_DB_t::XS &xs = d_mat->get(matid, columns[n]);

        // calculate the cdf for scattering to this group
        if (xs.has_upscatter(g))
            cdf += xs.sigma_s(g, 0) * total;

        // check to see if this is the group that the particle scatters to
        if (rnd <= cdf)
            return columns[n];
    }

    // now do downscatter groups
    int gp = 0;
    for (gp = g; gp < d_Ng; ++gp)
    {
        // calculate the cdf for scattering to this group
        cdf += d_mat->get(matid, gp).sigma_s(g, 0) * total;

        // see if we have sampled this group
        if (rnd <= cdf)
            return gp;
    }
    Check (nemesis::soft_equiv(cdf, 1.0));

    // we failed to sample
    Validate(false, "Failed to sample group.");
    return -1;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a fission group.
 *
 * This function is optimized based on the assumption that nearly all of the
 * fission emission is in the first couple of groups.
 */
template<class Geometry>
int MG_Physics<Geometry>::sample_fission_group(int    matid,
                                               double rnd) const
{
    Require (d_mat->assigned_fission(matid));

    // running cdf; we make the cdf on the fly because nearly all of the
    // emission is in the first couple of groups so its not worth storing for
    // a binary search
    double cdf = 0.0;

    // sample cdf
    for (int g = 0; g < d_Ng; ++g)
    {
        // update cdf
        cdf += d_mat->get_fission(matid).data(g, XS_DB_t::FIS_XS::CHI);

        // check for sampling; update particle's physics state and return
        if (rnd <= cdf)
        {
            // update the group in the particle
            return g;
        }
    }

    // we failed to sample
    Validate(false, "Failed to sample fission group.");
    return -1;
}

} // end namespace shift

#endif // mc_Physics_cc

//---------------------------------------------------------------------------//
//                 end of Physics.cc
//---------------------------------------------------------------------------//
