//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Physics.hh
 * \author Stuart Slattery
 * \brief  MG_Physics class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Physics_hh
#define cuda_mc_Physics_hh

#include <vector>
#include <memory>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "utils/Static_Map.hh"
#include "cuda_utils/Shared_Device_Ptr.hh"
#include "cuda_utils/CudaMacros.hh"
#include "cuda_utils/Stream.hh"
#include "xs/XS.hh"
#include "cuda_xs/XS_Device.hh"
#include "Definitions.hh"
#include "Group_Bounds.hh"
#include "Particle_Vector.hh"

namespace cuda_profugus
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
 * \example mc_physics/test/tstPhysics.cc
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

    typedef profugus::XS                        XS_t;
    typedef Teuchos::ParameterList              ParameterList_t;
    typedef typename Geometry_t::Space_Vector   Space_Vector;
    //@}

    //! Fission site structure for storing fission sites in k-code.
    struct Fission_Site
    {
        int          m;
        Space_Vector r;
    };

    //! Fission_Site container.
    typedef std::vector<Fission_Site> Fission_Site_Container;

    //! Fission site structure for storing fission sites in k-code on device.
    struct Device_Fission_Site
    {
        int          m;
        Space_Vector r;
	int          n;
    };
    
  private:
    // >>> DATA

    // Cross section database.
    cuda_utils::Shared_Device_Ptr<XS_Device> d_mat;

    // XS raw device pointer for inline device functions.
    XS_Device* d_mat_device;

    // Geometry.
    cuda_utils::Shared_Device_Ptr<Geometry> d_geometry;

  public:
    // Constructor that auto-creates group bounds.
    Physics( ParameterList_t& db, const XS_t& mat );

    // Destructor.
    ~Physics();

    // >>> PUBLIC TRANSPORT INTERFACE

    //! Set the geometry.
    void set_geometry(const cuda_utils::Shared_Device_Ptr<Geometry>& g)
    { REQUIRE(g); d_geometry = g; }

    //! Get the geometry.
    cuda_utils::Shared_Device_Ptr<Geometry> get_geometry() const { return d_geometry; }

    // Initialize the physics state.
    void initialize( const std::vector<double>& energy, 
		     cuda_utils::Shared_Device_Ptr<Particle_Vector_t>& particles );

    // Get a total cross section from the physics library.
    PROFUGUS_DEVICE_FUNCTION
    double total( const physics::Reaction_Type type,
		  const int matid,
		  const int group ) const
    {
	switch (type)
	{
	    case physics::TOTAL:
		return d_mat_device->vector(matid, XS_t::TOTAL)(group);
	    case physics::SCATTERING:
		return d_scatter[d_matid_g2l[matid] * d_mat_device->num_groups() + group];
	    case physics::FISSION:
		return d_mat_device->vector(matid, XS_t::SIG_F)(group);
	    case physics::NU_FISSION:
		return d_mat_device->vector(matid, XS_t::NU_SIG_F)(group);
	    default:
		return 0.0;
	}
    }

    //! Get the minimum energy allowed for a particle
    PROFUGUS_DEVICE_FUNCTION
    double min_energy() const { return d_gb_device->group_bounds()[d_Ng]; }

    //! Get the maximum energy allowed for a particle
    PROFUGUS_DEVICE_FUNCTION
    double max_energy() const { return d_gb_device->group_bounds()[0]; }

    // >>> TYPE-CONCEPT INTERFACE

    // Process particles through a collision.
    void collide( cuda_utils::Shared_Device_Ptr<Particle_Vector_t>& particles,
                  cuda_utils::Stream<cuda_utils::arch::Device> stream =
                  cuda_utils::Stream<cuda_utils::arch::Device>() );

    // Sample fission spectrum and initialize the physics state.
    PROFUGUS_DEVICE_FUNCTION
    void initialize_fission(const int matid,
			    const double ran,
			    int& group,
			    bool& sampled ) const
    {
	sampled = false;
	if (is_fissionable(matid))
	{
	    group = sample_fission_group( matid, ran );
	    sampled   = true;
	}
    }

    // Initialize a physics state at a fission site.
    PROFUGUS_DEVICE_FUNCTION
    void initialize_fission(const Fission_Site &fs,
			    const double ran,
			    int& group,
			    bool& sampled) const
    {
	REQUIRE( is_fissionable(fs.m) );
	group = sample_fission_group( fs.m, ran );
	sampled = true;
    }

    // Sample fission site.
    void sample_fission_site( cuda_utils::Shared_Device_Ptr<Particle_Vector_t>& particles,
                              Fission_Site_Container &fsc,
                              double keff,
                              cuda_utils::Stream<cuda_utils::arch::Device> stream =
                              cuda_utils::Stream<cuda_utils::arch::Device>() );

    // Return whether a given material is fissionable
    PROFUGUS_DEVICE_FUNCTION
    bool is_fissionable(unsigned int matid) const
    {
        return d_fissionable[d_matid_g2l[matid]];
    }

    // >>> FISSION SITE CONTAINER OPERATIONS

    //! Fission site position.
    static Space_Vector fission_site(const Fission_Site &fs) { return fs.r; }

    //! Size of fission sites in bytes.
    static int fission_site_bytes() { return sizeof(Fission_Site); }

    // >>> CLASS FUNCTIONS

    //! Get cross section database.
    const cuda_utils::Shared_Device_Ptr<XS_Device>& xs() const { return d_mat; }

    //! Number of discrete energy groups
    int num_groups() const { return d_Ng; }

    //! Group boundaries
    const cuda_utils::Shared_Device_Ptr<Group_Bounds>& group_bounds() const { return d_gb; }

  private:
    // >>> IMPLEMENTATION

    // Private types.
    typedef def::Vec_Dbl         Vec_Dbl;
    typedef def::Vec_Int         Vec_Int;
    typedef std::vector<Vec_Dbl> Vec_Vec_Dbl;

    // Boolean for implicit capture.
    bool d_implicit_capture;

    // Number of groups and materials in the material database.
    int d_Ng, d_Nm;

    // Group boundaries.
    cuda_utils::Shared_Device_Ptr<Group_Bounds> d_gb;

    // Raw pointer to group boundaries for inline device functions.
    Group_Bounds* d_gb_device;

    // Global-to-local mapping of matids. On-device.
    int* d_matid_g2l;

    // Total scattering for each group and material. On-device.
    double* d_scatter;

    // Fissionable bool by local matid. On-device.
    int* d_fissionable;

    // Host fission site work vector.
    std::vector<Device_Fission_Site> d_host_sites;
    
    // Device fission sites work vector.
    Device_Fission_Site* d_device_sites;

  private:

    //! Sample a fission group.
    PROFUGUS_DEVICE_FUNCTION
    int sample_fission_group(const unsigned int matid,
			     const double       rnd) const
    {
	double cdf = 0.0;
	const auto chi = d_mat_device->vector(matid, XS_t::CHI);
	for (int g = 0; g < d_mat_device->num_groups(); ++g)
	{
	    cdf += chi(g);
	    if (rnd <= cdf) return g;
	}
	return -1;
    }
};

} // end namespace cuda_profugus

#endif // cuda_mc_Physics_hh

//---------------------------------------------------------------------------//
//              end of cuda_mc/Physics.hh
//---------------------------------------------------------------------------//
