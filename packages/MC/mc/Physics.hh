//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Physics.hh
 * \author Thomas M. Evans
 * \date   Thursday May 1 11:15:7 2014
 * \brief  MG_Physics class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Physics_hh
#define mc_Physics_hh

#include <vector>
#include <memory>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "utils/Static_Map.hh"
#include "xs/XS.hh"
#include "geometry/Geometry.hh"
#include "Definitions.hh"
#include "Group_Bounds.hh"
#include "Particle.hh"
#include "Bank.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Physics
 * \brief Multigroup physics implementation for Monte Carlo transport.
 *
 * \section db_mg_physics Standard DB Entries for MG_Physics
 *
 * The following standard data entries control the MG physics:
 *
 * \arg \c implicit_capture (bool) does implicit capture (default: true)
 *
 * \arg \c neutron_bnd (vector<double>) neutron group boundaries
 *
 * \arg \c photon_bnd (vector<double>) photon group boundaries
 *
 * \arg \c check_balance (bool) check for balanced scattering tables (default:
 * false)
 *
 * Either the neutron boundaries, photon boundaries, or both \b must be defined.
 */
/*!
 * \example mc_physics/test/tstPhysics.cc
 *
 * Test of Physics.
 */
//===========================================================================//

class Physics
{
  public:
    //@{
    //! Useful typedefs.
    typedef Core                        Geometry_t;
    typedef Geometry_t::Geo_State_t     Geo_State_t;
    typedef Particle                    Particle_t;
    typedef Bank                        Bank_t;
    typedef std::shared_ptr<Geometry_t> SP_Geometry;
    typedef std::shared_ptr<Particle_t> SP_Particle;

    typedef Particle_t::RNG_t             RNG;
    typedef XS                            XS_t;
    typedef Teuchos::RCP<XS_t>            RCP_XS;
    typedef Teuchos::ParameterList        ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t> RCP_Std_DB;
    typedef Geometry_t::Space_Vector      Space_Vector;
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
    RCP_XS d_mat;

    // Geometry.
    SP_Geometry d_geometry;

  public:
    // Constructor that auto-creates group bounds.
    explicit Physics(RCP_Std_DB db, RCP_XS mat);

    // >>> PUBLIC TRANSPORT INTERFACE

    //! Set the geometry.
    void set_geometry(SP_Geometry g) { Require(g); d_geometry = g; }

    //! Get the geometry.
    SP_Geometry get_geometry() const { return d_geometry; }

    // Initialize the physics state.
    void initialize(double E, Particle_t &p);

    // Get a total cross section from the physics library.
    double total(physics::Reaction_Type type, unsigned int matid,
                 const Particle_t &p);

    //! Get the energy from a particle via its physics state
    double energy(const Particle_t &p) const
    {
        double low = 0.0, up = 0.0;
        d_gb.get_energy(p.group(), low, up);
        return low;
    }

    //! Get the minimum energy allowed for a particle
    double min_energy() const { return d_gb.group_bounds()[d_Ng]; }

    //! Get the maximum energy allowed for a particle
    double max_energy() const { return d_gb.group_bounds()[0]; }

    // >>> TYPE-CONCEPT INTERFACE

    // Process a particle through a physical collision.
    void collide(Particle_t &particle, Bank_t &bank);

    // Sample fission site.
    int sample_fission_site(const Particle_t &p, Fission_Site_Container &fsc,
                            double keff);

    // Sample fission spectrum and initialize the physics state.
    bool initialize_fission(unsigned int matid, Particle_t &p);

    // Initialize a physics state at a fission site.
    bool initialize_fission(Fission_Site &fs, Particle_t &p);

    // Return whether a given material is fissionable
    bool is_fissionable(unsigned int matid) const
    {
        return d_fissionable[d_mid2l[matid]];
    }

    // >>> FISSION SITE CONTAINER OPERATIONS

    //! Fission site position.
    static Space_Vector fission_site(const Fission_Site &fs) { return fs.r; }

    //! Size of fission sites in bytes.
    static int fission_site_bytes() { return sizeof(Fission_Site); }

    // >>> CLASS FUNCTIONS

    //! Get cross section database.
    const XS_t& xs() const { return *d_mat; }

    //! Number of discrete energy groups
    int num_groups() const { return d_Ng; }

    //! Group boundaries
    const Group_Bounds& group_bounds() const { return d_gb; }

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
    Group_Bounds d_gb;

    // Matid-to-local hash such that d_mid2l[matid] = [0,N).
    Static_Map<unsigned int, unsigned int> d_mid2l;

    // Total scattering for each group and material.
    Vec_Vec_Dbl d_scatter;

    // Fissionable bool by local matid.
    std::vector<bool> d_fissionable;

    // Material id of current region.
    int d_matid;

    // Sample a group.
    int sample_group(int matid, int g, double rnd) const;

    // Sample a fission group.
    int sample_fission_group(unsigned int matid, double rnd) const;
};

} // end namespace shift

#endif // mc_Physics_hh

//---------------------------------------------------------------------------//
//              end of mc/Physics.hh
//---------------------------------------------------------------------------//
