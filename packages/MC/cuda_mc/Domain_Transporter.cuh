//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Domain_Transporter.cuh
 * \author Thomas M. Evans
 * \date   Mon May 12 12:02:13 2014
 * \brief  Domain_Transporter class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Domain_Transporter_cuh
#define cuda_mc_Domain_Transporter_cuh

#include <memory>
#include <thrust/device_vector.h>

#include "Particle_Vector.cuh"
#include "Physics.cuh"
#include "Step_Selector.cuh"
#include "VR_Roulette.cuh"
#include "Tallier.cuh"
#include "Definitions.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Domain_Transporter
 * \brief Transport a particle on a computational domain.
 *
 * This class does no communication; it takes a particle and transports it
 * until it leaves the domain.
 */
/*!
 * \example mc/test/tstDomain_Transporter.cc
 *
 * Test of Domain_Transporter.
 */
//===========================================================================//

template <class Geometry>
class Domain_Transporter
{
  public:
    //@{
    //! Useful typedefs.
    typedef Geometry                                   Geometry_t;
    typedef Physics<Geometry_t>                        Physics_t;
    typedef cuda_utils::Space_Vector                   Space_Vector;
    typedef Particle_Vector<Geometry_t>                Particle_Vector_t;
    typedef VR_Roulette<Geometry_t>                    VR_Roulette_t;
    typedef Tallier<Geometry_t>                        Tallier_t;
    typedef Teuchos::RCP<Teuchos::ParameterList>       RCP_Std_DB;
    typedef thrust::device_vector<Fission_Site>        Fission_Site_Vector;
    typedef std::shared_ptr<Fission_Site_Vector>       SP_Fission_Site_Vec;
    //@}

  private:
    // >>> DATA

    // Problem geometry implementation.
    Geometry    *d_geometry;

    // Problem physics implementation.
    Physics_t  *d_physics;

    // Regular tallies.
    Tallier_t *d_tallier;

    // Variance reduction.
    VR_Roulette_t *d_vr;

    // Fission sites.
    Fission_Site *d_fission_sites;

  public:

    // Constructor.
    Domain_Transporter(Geometry *geometry,
                       Physics_t *physics,
                       Tallier_t *tallier,
                       VR_Roulette_t *vr,
                       Fission_Site  *fission_sites,
                       bool           sample_fission_sites,
                       int           *num_fission_sites,
                       int            max_fission_sites,
                       double         keff,
                       int            max_steps)
        : d_geometry(geometry)
        , d_physics(physics)
        , d_tallier(tallier)
        , d_vr(vr)
        , d_fission_sites(fission_sites)
        , d_sample_fission_sites(sample_fission_sites)
        , d_num_fission_sites(num_fission_sites)
        , d_max_fission_sites(max_fission_sites)
        , d_keff(keff)
        , d_max_steps(max_steps)
    {
    }

    // Transport a particle through the domain.
    __device__ void transport(int pid, Particle_Vector_t &particles) const;

    //! Return the number of sampled fission sites.
    __device__
    int num_sampled_fission_sites() const { return *d_num_fission_sites; }

  private:
    // >>> IMPLEMENTATION

    // Flag indicating that fission sites should be sampled.
    bool d_sample_fission_sites;

    // Number of fission sites sampled.
    int *d_num_fission_sites;

    // Number of fission sites allocated
    int d_max_fission_sites;

    // Current keff iterate.
    double d_keff;

    // Number of steps per transport solve
    int d_max_steps;

    // Process collisions and boundaries.
    __device__ void process_boundary(int pid, Particle_Vector_t &particles) const;
    __device__ void process_collision(int pid, Particle_Vector_t &particles,
                                      double      step) const;
};

//===========================================================================//
/*!
 * \class Domain_Transporter_DMM
 * \brief Device memory manager for Domain_Transporter.
 */
//===========================================================================//

template <class Geometry>
class Domain_Transporter_DMM :
    public cuda::Device_Memory_Manager<Domain_Transporter<Geometry>>
{
  public:
    //@{
    //! Useful typedefs.
    typedef Geometry                                   Geometry_t;
    typedef Physics<Geometry_t>                        Physics_t;
    typedef cuda_utils::Space_Vector                   Space_Vector;
    typedef VR_Roulette<Geometry_t>                    VR_Roulette_t;
    typedef Tallier<Geometry_t>                        Tallier_t;
    typedef Teuchos::RCP<Teuchos::ParameterList>       RCP_Std_DB;
    typedef thrust::device_vector<Fission_Site>        Fission_Site_Vector;
    typedef std::shared_ptr<Fission_Site_Vector>       SP_Fission_Site_Vec;
    //@}

    //@{
    //! Smart pointers.
    typedef cuda::Shared_Device_Ptr<Geometry_t>       SDP_Geometry;
    typedef cuda::Shared_Device_Ptr<Physics_t>        SDP_Physics;
    typedef cuda::Shared_Device_Ptr<VR_Roulette_t>    SDP_VR;
    typedef cuda::Shared_Device_Ptr<Tallier_t>        SDP_Tallier;
    //@}

  private:
    // >>> DATA

    // Problem geometry implementation.
    SDP_Geometry d_geometry;

    // Problem physics implementation.
    SDP_Physics d_physics;

    // Regular tallies.
    SDP_Tallier d_tallier;

    // Variance reduction.
    SDP_VR d_vr;

    // Fission sites.
    SP_Fission_Site_Vec d_fission_sites;

  public:

    // Constructor.
    Domain_Transporter_DMM(RCP_Std_DB   db,
                           SDP_Geometry geometry,
                           SDP_Physics  physics,
                           SDP_VR       vr      = SDP_VR() );

    // DMM interface
    Domain_Transporter<Geometry> device_instance()
    {
        // Get pointer to fission sites
        // Assign null ptr if not allocated
        Fission_Site *sites = nullptr;
        if (d_fission_sites)
        {
            REQUIRE(d_sample_fission_sites);
            sites = d_fission_sites->data().get();
        }
        else
        {
            REQUIRE(!d_sample_fission_sites);
            REQUIRE(d_num_fission_sites.empty());
        }

        Domain_Transporter<Geometry> transporter(
                d_geometry.get_device_ptr(),
                d_physics.get_device_ptr(),
                d_tallier.get_device_ptr(),
                d_vr.get_device_ptr(),
                sites,
                d_sample_fission_sites,
                d_num_fission_sites.data().get(),
                d_max_fission_sites,
                d_keff,
                d_max_steps);
        return transporter;
    }

    // Set tallier
    void set(SDP_Tallier tallier);

    // Set fission site sampling.
    void set(SP_Fission_Site_Vec fission_sites, double keff);

    //! Return the number of sampled fission sites.
    int num_sampled_fission_sites() const { return d_num_fission_sites[0]; }

  private:
    // >>> IMPLEMENTATION

    // Total cross section in region.
    double d_xs_tot;

    // Flag indicating that fission sites should be sampled.
    bool d_sample_fission_sites;

    // Number of fission sites sampled.
    thrust::device_vector<int> d_num_fission_sites;

    // Number of fission sites allocated
    int d_max_fission_sites;

    // Current keff iterate.
    double d_keff;

    // Number of steps per transport solve
    int d_max_steps;
};

} // end namespace cuda_mc 

#include "Domain_Transporter.i.cuh"

#endif // cuda_mc_Domain_Transporter_cuh


//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.cuh
//---------------------------------------------------------------------------//
