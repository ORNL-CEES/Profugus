//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_physics/MG_Physics.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 04 12:50:33 2011
 * \brief  MG_Physics class definition.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_physics_MG_Physics_hh
#define mc_physics_MG_Physics_hh

#include "Physics.hh"

#include <vector>

#include "utils/SP.hh"
#include "utils/Definitions.hh"
#include "material/XS_DB.hh"
#include "material/S_Graph.hh"
#include "database/Std_DB.hh"
#include "mc/Definitions.hh"
#include "mc/Group_Bounds.hh"

#include "Particle.hh"
#include "Shift_Bank.hh"
#include "Continuous_Fission_Comm.hh"
#include "Continuous_Fission_Rebalance.hh"
#include "MG_State.hh"

namespace shift
{

//===========================================================================//
/*!
 * \class MG_Physics
 * \brief Multigroup physics implementation for Monte Carlo transport.
 *
 * The derived physics class, MG_Physics, must define all of the functions
 * declared in shift::Physics.  The following additional types are required:
 * \code
   typedef MG_Physics<Geo_State>                      Physics_t;
   typedef typename Physics_t::Fission_Site           Fission_Site;
   typedef typename Physics_t::Fission_Site_Container Fission_Site_Container;
   typedef typename Physics_t::Fission_Comm_t         Fission_Comm_t;
   typedef typename Physics_t::SP_Fission_Comm        SP_Fission_Comm;
   typedef typename Physics_t::Fission_Rebalance_t    Fission_Rebalance_t;
   typedef typename Physics_t::SP_Fission_Rebalance   SP_Fission_Rebalance;
 * \endcode
 * The following functions are required as part of the Physics-type concept:
 * \code

   void Physics_t::sample_fission_site(const Particle_t       &p,
                                       Fission_Site_Container &fsc,
                                       double                  keff);

   bool Physics_t::initialize_fission(int matid, RNG rng,
                                      Physics_State_t &state);

   void Physics_t::initialize(const Fission_Site &fs,
                              Physics_State_t    &state);

   int num_fission_particles(const Fission_Site &fs) const;

   void subtract_fission_particle(Fission_Site &fs) const;

   Space_Vector fission_site(const Fission_Site &fs) const;
 * \endcode
 *
 * The requirements on the fission site container are defined by
 * Continuous_Fission_Comm and Continuous_Fission_Rebalance.
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
 * \example mc_physics/test/tstMG_Physics.cc
 *
 * Test of MG_Physics.
 */
//===========================================================================//

template<class Geometry>
class MG_Physics : public Physics<Geometry, MG_State>
{
    typedef Physics<Geometry, MG_State> Base;
    typedef MG_Physics<Geometry>        Physics_t;

  public:
    //@{
    //! Useful typedefs.
    typedef Geometry                               Geometry_t;
    typedef MG_State                               Physics_State_t;
    typedef typename Base::Geo_State_t             Geo_State_t;
    typedef Particle<Geo_State_t, Physics_State_t> Particle_t;
    typedef shift::Bank<Geometry_t, Physics_t>     Bank_t;
    typedef nemesis::SP<Geometry_t>                SP_Geometry;
    typedef nemesis::SP<Particle_t>                SP_Particle;

    typedef typename Particle_t::RNG               RNG;
    typedef denovo::XS_DB                          XS_DB_t;
    typedef nemesis::SP<XS_DB_t>                   SP_XS_DB;
    typedef nemesis::SP<database::Std_DB>          SP_Std_DB;
    typedef nemesis::SP<mc::Group_Bounds>          SP_Group_Bounds;
    typedef typename Geometry_t::Space_Vector      Space_Vector;
    //@}

    //! Flag for being multigroup in energy
    typedef mc::physics::MG Energy_Treatment_t;

    //! Fission site structure for storing fission sites in k-code.
    struct Fission_Site
    {
        int          m;
        Space_Vector r;
    };

    //! Fission_Site container.
    typedef std::vector<Fission_Site> Fission_Site_Container;

    //@{
    //! Fission site communicator.
    typedef Continuous_Fission_Comm< MG_Physics<Geometry> > Fission_Comm_t;
    typedef nemesis::SP<Fission_Comm_t>                     SP_Fission_Comm;
    //@}

    //@{
    //! Fission site rebalance.
    typedef Continuous_Fission_Rebalance<Physics_t> Fission_Rebalance_t;
    typedef nemesis::SP<Fission_Rebalance_t>        SP_Fission_Rebalance;
    //@}

  private:
    // >>> DATA

    // Cross section database.
    SP_XS_DB d_mat;

    // Make protected geometry available.
    using Base::b_geometry;

  public:
    // Constructor that auto-creates group bounds
    explicit MG_Physics(SP_Std_DB db, SP_XS_DB mat);

    // Constructor.
    explicit MG_Physics(SP_Std_DB db, SP_XS_DB mat, SP_Group_Bounds gb);

    // >>> DERIVED PUBLIC INTERFACE

    // No need to retain RNG outside of particle
    void set_rng(const RNG& rng) { /* * */ }

    // Initialize the physics state.
    void initialize(mc::Particle_Type type, double E, Physics_State_t &state);

    // Get a total cross section from the physics library.
    double total(mc::physics::Reaction_Type type, int matid,
                 const Physics_State_t &state);

    //! Get the energy from a particle via its physics state
    double energy(const Physics_State_t& state) const
    {
        Require(state.group < d_Ng);
        return d_gb->get_lower_energy(state.group);
    }

    //! Get the minimum energy allowed for a particle based on type
    double min_energy(mc::Particle_Type type) const;

    //! Get the maximum energy allowed for a particle based on type
    double max_energy(mc::Particle_Type type) const;

    //! Get the particle's type (photon, neutron) via its physics state
    mc::Particle_Type particle_type(const Physics_State_t &state) const
    {
        return state.type;
    }

    // >>> TYPE-CONCEPT INTERFACE

    // Process a particle through a physical collision.
    void collide(Particle_t &particle, Bank_t &bank);

    // Sample fission site.
    int sample_fission_site(const Particle_t &particle,
                            Fission_Site_Container &fsc, double keff);

    // Sample fission spectrum and initialize the physics state.
    bool initialize_fission(int matid, RNG rng, Physics_State_t &state);

    // Initialize a physics state at a fission site.
    bool initialize_fission(Fission_Site &fs, RNG rng, Physics_State_t &state);

    // Return whether a given material is fissionable
    bool is_fissionable(int matid) const
    {
        return d_mat->assigned_fission(matid) > 0 ? true : false;
    }

    //! Return false because MG physics does not inherently split particles
    bool uses_splitting() const { return false; }

    //! Pickle the physics state.
    void pickle(Physics_State_t &state) { /* * */ }

    //! Restore the physics from a persistent state.
    void restore(Physics_State_t &pickled_state) { /* * */}

    //! Called after transport solve; save diagnostics and clear
    void finalize() { /* * */ }

    //! Rebuild internal data before a new depletion step
    void reset() { /* * */ }

    // >>> FISSION SITE CONTAINER OPERATIONS

    //! Fission site position.
    static Space_Vector fission_site(const Fission_Site &fs) { return fs.r; }

    //! Size of fission sites in bytes.
    static int fission_site_bytes() { return sizeof(Fission_Site); }

    // >>> CLASS FUNCTIONS

    //! Get cross section database.
    const XS_DB_t& xsdb() const { return *d_mat; }

    //! Number of discrete energy groups
    int num_groups() const { return d_Ng; }

    //! Group boundaries
    SP_Group_Bounds group_bounds() const { return d_gb; }

  private:
    // >>> IMPLEMENTATION

    // Private types.
    typedef def::Vec_Dbl         Vec_Dbl;
    typedef def::Vec_Int         Vec_Int;
    typedef std::vector<Vec_Dbl> Vec_Vec_Dbl;
    typedef denovo::S_Graph      S_Graph;

    // Boolean for implicit capture.
    bool d_implicit_capture;

    // Check cross section balance.
    bool d_check_balance;

    // Number of groups and materials in the material database.
    int d_Ng, d_Nm;

    // Group boundaries.
    SP_Group_Bounds d_gb;

    // Total scattering for each group and material.
    Vec_Vec_Dbl d_scatter;

    // Graph of scattering matrix.
    S_Graph d_graph;

    // Material id of current region.
    int d_matid;

    // Construct private data
    void construct(SP_Std_DB db, SP_XS_DB mat, SP_Group_Bounds gb);

    // Sample a group.
    int sample_group(int matid, int g, double rnd) const;

    // Sample a fission group.
    int sample_fission_group(int matid, double rnd) const;
};

} // end namespace shift

#endif // mc_physics_MG_Physics_hh

//---------------------------------------------------------------------------//
//              end of mc_physics/MG_Physics.hh
//---------------------------------------------------------------------------//
