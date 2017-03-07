//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Physics.cuh
 * \author Thomas M. Evans
 * \date   Thursday May 1 11:15:7 2014
 * \brief  MG_Physics class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Physics_cuh
#define cuda_mc_Physics_cuh

#include <vector>
#include <memory>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CudaUtils/cuda_utils/Shared_Device_Ptr.hh"
#include "Matprop/cuda_xs/XS_Device.hh"
#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "xs/XS.hh"
#include "mc/Definitions.hh"
#include "Particle_Vector.cuh"
#include "Definitions.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Physics
 * \brief Multigroup physics implementation for Monte Carlo transport.
 *
 * \section db_physics Standard DB Entries for Physics
 *
 * The following standard data entries control the physics:
 *
 * \arg \c implicit_capture (bool) does implicit capture (default: true)
 *
 * \arg \c check_balance (bool) check for balanced scattering tables (default:
 * false)
 */
/*!
 * \example cuda_mc/test/tstPhysics.cc
 *
 * Test of Physics.
 */
//===========================================================================//

template <class Geometry>
class Physics
{
  public:
    //@{
    //! Useful typedefs.
    typedef Geometry                            Geometry_t;
    typedef typename Geometry_t::Geo_State_t    Geo_State_t;
    typedef Particle_Vector<Geometry_t>         Particle_Vector_t;
    typedef cuda::Shared_Device_Ptr<Geometry_t> SDP_Geometry;
    typedef cuda_profugus::XS_Device            XS_Dev_t;
    typedef cuda::Shared_Device_Ptr<XS_Dev_t>   SDP_XS_Dev;
    typedef profugus::XS                        XS_t;
    typedef Teuchos::RCP<XS_t>                  RCP_XS;
    typedef Teuchos::ParameterList              ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t>       RCP_Std_DB;
    typedef cuda_utils::Space_Vector            Space_Vector;
    typedef profugus::physics::Reaction_Type    Reaction_Type;
    //@}

  private:
    // >>> DATA

    // Cross section database.
    XS_Dev_t  *d_mat;

    // Geometry.
    Geometry    *d_geometry;

  public:
    // Constructor that auto-creates group bounds.
    explicit Physics(RCP_Std_DB db, RCP_XS mat_host, SDP_XS_Dev mat);

    // Due to Cuda memory management, disallow copy and assignment
    Physics(const Physics &phys)            = delete;
    Physics& operator=(const Physics &phys) = delete;

    ~Physics();

    // >>> PUBLIC TRANSPORT INTERFACE

    //! Set the geometry.
    void set_geometry(SDP_Geometry g)
    {
        REQUIRE(g.get_host_ptr());
        REQUIRE(g.get_device_ptr());
        d_geometry = g.get_device_ptr();
    }

    // Get a total cross section from the physics library.
    __device__ double total(Reaction_Type type, int pid,
                            const Particle_Vector_t &particles) const;

    // >>> TYPE-CONCEPT INTERFACE

    // Process a particle through a physical collision.
    __device__ void collide(int pid, Particle_Vector_t &particles) const;

    // Sample fission site.
    __device__ int sample_fission_site(int                        pid,
                                       const Particle_Vector_t   &particles,
                                             double               keff) const;

    // Sample fission spectrum and initialize the physics state.
    __device__ bool initialize_fission(unsigned int matid, int pid,
                                       Particle_Vector_t &particles) const;

    // Initialize a physics state at a fission site.
    __device__ bool initialize_fission(const Fission_Site &fs,
                                       int                 pid,
                                       Particle_Vector_t  &particles) const;

    // Return whether a given material is fissionable
    __device__ bool is_fissionable(unsigned int matid) const
    {
        return d_fissionable[matid];
    }

    // >>> FISSION SITE CONTAINER OPERATIONS

    // >>> CLASS FUNCTIONS

    //! Get cross section database (device).
    __device__ XS_Dev_t  * xs_device() const { return d_mat; }

    //! Number of discrete energy groups
    __host__ __device__ int num_groups() const { return d_Ng; }

  private:
    // >>> IMPLEMENTATION

    // Boolean for implicit capture.
    bool d_implicit_capture;

    // Number of groups and materials in the material database.
    int d_Ng, d_Nm;

    // Total scattering for each group and material.
    double *d_scatter;

    // Fissionable bool by local matid (stored as int for host/device
    // compatibility)
    int *d_fissionable;

    // Index for a material/group combination
    __host__ __device__ int group_mat_index(int g, int m) const
    {
        DEVICE_REQUIRE( g < d_Ng );
        DEVICE_REQUIRE( m < d_Nm );
        return g + d_Ng * m;
    }

    // Sample a group.
    __device__ int sample_group(int matid, int g, double rnd) const;

    // Sample a fission group.
    __device__ int sample_fission_group(unsigned int matid, double rnd) const;
};

} // end namespace cuda_mc

#include "Physics.i.cuh"

#endif // cuda_mc_Physics_cuh

//---------------------------------------------------------------------------//
//              end of cuda_mc/Physics.cuh
//---------------------------------------------------------------------------//
