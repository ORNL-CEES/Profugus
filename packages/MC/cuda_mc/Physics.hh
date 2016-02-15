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

    typedef XS                                  XS_t;
    typedef Teuchos::RCP<XS_t>                  RCP_XS;
    typedef Teuchos::ParameterList              ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t>       RCP_Std_DB;
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

  private:
    // >>> DATA

    // Cross section database.
    Shared_Device_Ptr<XS_Device> d_mat;

    // Geometry.
    Shared_Device_Ptr<Geometry> d_geometry;

  public:
    // Constructor that auto-creates group bounds.
    explicit Physics(RCP_Std_DB db, RCP_XS mat);

    // Destructor.
    ~Physics();

    // >>> PUBLIC TRANSPORT INTERFACE

    //! Set the geometry.
    void set_geometry(const Shared_Device_Ptr<Geometry>& g) { REQUIRE(g); d_geometry = g; }

    //! Get the geometry.
    Shared_Device_Ptr<Geometry> get_geometry() const { return d_geometry; }

    // Initialize the physics state.
    void initialize(const std::vector<double>& energy, 
		    Shared_Device_Ptr<Particle_Vector_t>& particles );

    // Get a total cross section from the physics library.
    double total(physics::Reaction_Type type, const Particle_t &p);

    //! Get the energy from a particle via its physics state
    double energy(const Particle_t &p) const;

    //! Get the minimum energy allowed for a particle
    PROFUGUS_DEVICE_FUNCTION
    double min_energy() const { return d_gb.group_bounds()[d_Ng]; }

    //! Get the maximum energy allowed for a particle
    PROFUGUS_DEVICE_FUNCTION
    double max_energy() const { return d_gb.group_bounds()[0]; }

    // >>> TYPE-CONCEPT INTERFACE

    // Process a particle through a physical collision.
    void collide(Particle_t &particle);

    // Sample fission site.
    int sample_fission_site(const Particle_t &p, Fission_Site_Container &fsc,
                            double keff);

    // Sample fission spectrum and initialize the physics state.
    bool initialize_fission(unsigned int matid, Particle_t &p);

    // Initialize a physics state at a fission site.
    bool initialize_fission(Fission_Site &fs, Particle_t &p);

    // Return whether a given material is fissionable
    PROFUGUS_DEVICE_FUNCTION
    bool is_fissionable(unsigned int matid) const
    {
        return d_fissionable[d_mid_g2l[matid]];
    }

    // >>> FISSION SITE CONTAINER OPERATIONS

    //! Fission site position.
    static Space_Vector fission_site(const Fission_Site &fs) { return fs.r; }

    //! Size of fission sites in bytes.
    static int fission_site_bytes() { return sizeof(Fission_Site); }

    // >>> CLASS FUNCTIONS

    //! Get cross section database.
    const Shared_Device_Ptr<XS>& xs() const { return *d_mat; }

    //! Number of discrete energy groups
    int num_groups() const { return d_Ng; }

    //! Group boundaries
    const Shared_Device_Ptr<Group_Bounds>& group_bounds() const { return d_gb; }

  private:
    // >>> IMPLEMENTATION

    // Private types.
    typedef def::Vec_Dbl         Vec_Dbl;
    typedef def::Vec_Int         Vec_Int;
    typedef std::vector<Vec_Dbl> Vec_Vec_Dbl;

    // Boolean for implicit capture.
    bool d_implicit_capture;

    // Check cross section balance.
    bool d_check_balance;

    // Number of groups and materials in the material database.
    int d_Ng, d_Nm;

    // Group boundaries.
    Shared_Device_Ptr<Group_Bounds> d_gb;

    // Global-to-local mapping of matids. On-device.
    int* d_matid_g2l;

    // Total scattering for each group and material. On-device.
    double* d_scatter;

    // Fissionable bool by local matid. On-device.
    int* d_fissionable;

    // Material id of current region.
    int d_matid;

    // Sample a group.
    int sample_group(int matid, int g, double rnd) const;

    // Sample a fission group.
    int sample_fission_group(unsigned int matid, double rnd) const;
};

} // end namespace cuda_profugus

#endif // cuda_mc_Physics_hh

//---------------------------------------------------------------------------//
//              end of cuda_mc/Physics.hh
//---------------------------------------------------------------------------//
