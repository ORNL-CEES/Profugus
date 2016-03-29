//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/KDE_Fission_Source.cc
 * \author Gregory G. Davidson
 * \date   Mon Nov 23 15:47:23 2015
 * \brief  KDE_Fission_Source class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "KDE_Fission_Source.hh"

#include <cmath>

#include "Axial_KDE_Kernel.hh"
#include "comm/Timing.hh"
#include "mc/Global_RNG.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
KDE_Fission_Source::KDE_Fission_Source(RCP_Std_DB     db,
                                       SP_Geometry    geometry,
                                       SP_Physics     physics,
                                       SP_RNG_Control rng_control)
    : Base(db, geometry, physics, rng_control)
{
    REQUIRE(db.is_valid_ptr());
    REQUIRE(db->isSublist("kde_db"));

    // Get the KDE database
    auto &kde_db = db->sublist("kde_db");

    VALIDATE(kde_db.isType<std::string>("kernel_type"),
             "KDE Kernel type was not specified in input");

    // Get the type of KDE kernel
    const std::string &kernel_type = kde_db.get<std::string>("kernel_type");

    // Get the bandwidth coefficient
    double coeff    = kde_db.get<double>("bnd_coeff", 1.06);
    double exponent = kde_db.get<double>("bnd_exp", -0.20);

    if (kernel_type == "fission_rejection")
    {
        // Instantiate the resample kernel
        d_kernel = std::make_shared<Axial_KDE_Kernel>(
            geometry, physics, Axial_KDE_Kernel::FISSION_REJECTION,
            coeff, exponent);
    }
    else if (kernel_type == "cell_rejection")
    {
        // Instantiate the resample kernel
        d_kernel = std::make_shared<Axial_KDE_Kernel>(
            geometry, physics, Axial_KDE_Kernel::CELL_REJECTION,
            coeff, exponent);
    }
    else
    {
        VALIDATE(false, "Unrecognized kernel type " << kernel_type);
    }

    ENSURE(d_kernel);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build a source from a fission site container.
 */
void KDE_Fission_Source::build_source(SP_Fission_Sites &fission_sites)
{
    REQUIRE(d_kernel);

    // Call Base method first.  Afterwards, we must use the internal fission
    // site container
    Base::build_source(fission_sites);
    CHECK(fission_sites->empty());

    // Calculate the bandwidths
    d_kernel->calc_bandwidths(*d_fission_sites);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a particle.
 */
KDE_Fission_Source::SP_Particle
KDE_Fission_Source::get_particle()
{
    using def::I; using def::J; using def::K;

    REQUIRE(d_wt > 0.0);
    REQUIRE(profugus::Global_RNG::d_rng.assigned());

    // particle
    SP_Particle p;
    CHECK(!p);

    // return a null particle if no source
    if (!d_num_left)
    {
        ENSURE(d_fission_sites ? d_fission_sites->empty() : true);
        return p;
    }

    SCOPED_TIMER_2("MC::KDE_Fission_Source.get_particle");

    // make a particle
    p = std::make_shared<Particle_t>();

    // use the global rng on this domain for the random number generator
    p->set_rng(profugus::Global_RNG::d_rng);
    RNG rng = p->rng();

    // material id
    int matid = 0;

    // particle position and isotropic direction
    Space_Vector r, omega;

    // sample the angle isotropically
    Base::sample_angle(omega, rng);

    // sample flag
    bool sampled;

    // if there is a fission site container than get the particle from there;
    // otherwise assume this is an initial source
    if (!is_initial_source())
    {
        CHECK(!d_fission_sites->empty());

        // Sample the fission site randomly
        Fission_Site &fs = this->sample_fission_site(rng);

        // get the location of the physics site
        r = b_physics->fission_site(fs);

        // Now, sample from the kernel
        r = d_kernel->sample_position(r, rng);

        // intialize the geometry state
        b_geometry->initialize(r, omega, p->geo_state());

        // get the material id
        matid = b_geometry->matid(p->geo_state());

        // initialize the physics state at the fission site
        sampled = b_physics->initialize_fission(fs, *p);
        CHECK(sampled);
    }
    else
    {
        matid = sample_geometry(r, omega, *p, rng);
    }

    // set the material id in the particle
    p->set_matid(matid);

    // set particle weight
    p->set_wt(d_wt);

    // make particle alive
    p->live();

    // update counters
    d_num_left--;
    d_num_run++;

    ENSURE(p->matid() == matid);
    return p;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a fission site from the fission-site population.
 */
KDE_Fission_Source::Fission_Site&
KDE_Fission_Source::sample_fission_site(RNG &rng) const
{
    // Get the number of fission sites
    unsigned int num_fs = d_fission_sites->size();

    // Sample a number between 0 and num_fs and floor it.  Reject in the
    // unlikely event that random number is 1.0
    unsigned int index;
    do
    {
        index = static_cast<unsigned int>(std::floor(rng.ran()*num_fs));
    } while(index == num_fs);
    CHECK(index < num_fs);

    return (*d_fission_sites)[index];
}

//---------------------------------------------------------------------------//
} // end namespace profugus

//---------------------------------------------------------------------------//
// end of MC/mc/KDE_Fission_Source.cc
//---------------------------------------------------------------------------//
