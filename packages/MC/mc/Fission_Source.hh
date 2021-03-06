//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/Fission_Source.hh
 * \author Thomas M. Evans
 * \date   Mon May 05 14:22:46 2014
 * \brief  Fission_Source class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Fission_Source_hh
#define MC_mc_Fission_Source_hh

#include "Teuchos_ArrayView.hpp"

#include "geometry/Cartesian_Mesh.hh"
#include "Fission_Rebalance.hh"
#include "Source.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Fission_Source
 * \brief Defines both an initial and within-cycle fission source for k-code
 * calculations.
 *
 * Through the database passed into the fission source the client specifies
 * the number of fission particles to run each cycle.  On cycles other than
 * the first, the number of fission source particles is determined by the
 * number of sites sampled during the previous generation.  To keep the
 * effective number of particles constant between cycles, fission particles
 * are born with the following weight,
 * \f[
   w = N_p / M\:,
 * \f]
 * where \f$N_p\f$ is the number of requested particles per cycle and \f$M\f$
 * is the total number of particles sampled from the fission source.  In the
 * first cycle \f$M = N_p\f$; however, in various parallel decompositions the
 * number of particles run in the first cycle may be slightly less than the
 * number of particles requested (but never more).  Thus, even in the first
 * cycle the weight of particles may be slightly greater than one because the
 * weight per particle is normalized to the \b requested number of particles
 * per cycle.  Using this weight preserves the constant total weight per
 * cycle.  The requested number of particles can be queried through the Np()
 * member function.
 *
 * \section fission_source_db Standard DB Entries for Fission_Source
 *
 * The database entries that control the profugus::Fission_Source are:
 *
 * \arg \c init_fission_src (vector<double>) box defined as (lox, hix, loy,
 * hiy, loz, hiz) in which the initial fission source is sampled (default: the
 * entire geometry)
 *
 * \arg \c Np (int) number of particles to use in each cycle (default:
 * 1000)
 */
/*!
 * \example mc/test/tstFission_Source.cc
 *
 * Test of Fission_Source.
 */
//===========================================================================//

template <class Geometry>
class Fission_Source : public Source<Geometry>
{
  public:
    //@{
    //! Typedefs.
    typedef Geometry                                    Geometry_t;
    typedef Physics<Geometry_t>                         Physics_t;
    typedef typename Physics_t::RCP_Std_DB              RCP_Std_DB;
    typedef typename Physics_t::Fission_Site            Fission_Site;
    typedef typename Physics_t::Fission_Site_Container  Fission_Site_Container;
    typedef typename Physics_t::Particle_t              Particle_t;
    typedef typename Geometry_t::Space_Vector           Space_Vector;
    typedef RNG_Control                                 RNG_Control_t;
    typedef RNG_Control_t::RNG_t                        RNG_t;
    typedef std::shared_ptr<Geometry>                   SP_Geometry;
    typedef std::shared_ptr<Physics_t>                  SP_Physics;
    typedef std::shared_ptr<Particle_t>                 SP_Particle;
    typedef std::shared_ptr<RNG_Control>                SP_RNG_Control;
    typedef std::shared_ptr<Fission_Site_Container>     SP_Fission_Sites;
    typedef Fission_Rebalance<Geometry_t>               Fission_Rebalance_t;
    typedef std::shared_ptr<Fission_Rebalance_t>        SP_Fission_Rebalance;
    typedef std::shared_ptr<Cartesian_Mesh>             SP_Cart_Mesh;
    typedef Teuchos::ArrayView<const double>            Const_Array_View;
    typedef def::Vec_Dbl                                Vec_Dbl;
    typedef def::Vec_Int                                Vec_Int;
    typedef def::size_type                              size_type;
    //@}

  protected:
    // >>> DATA

    // Fission site container.
    SP_Fission_Sites d_fission_sites;

    // Fission rebalance (across sets).
    SP_Fission_Rebalance d_fission_rebalance;

  public:
    // Constructor.
    Fission_Source(RCP_Std_DB db, SP_Geometry geometry, SP_Physics physics,
                   SP_RNG_Control rng_control);

    // Build the initial fission source.
    void build_initial_source();

    // Build the initial source from a mesh distribution.
    void build_initial_source(SP_Cart_Mesh mesh, Const_Array_View fis_dens);

    // Build a source from a fission site container.
    virtual void build_source(SP_Fission_Sites &fission_sites);

    // Create a fission site container.
    SP_Fission_Sites create_fission_site_container() const;

    // >>> DERIVED PUBLIC INTERFACE

    // Get a particle from the source.
    virtual SP_Particle get_particle();

    //! Boolean operator for source (true when source still has particles).
    bool empty() const { return d_num_left == 0; }

    //! Whether this is the initial source distribution or is unbuilt
    bool is_initial_source() const { return !d_fission_sites; }

    //! Number of particles to transport in the source on the current domain.
    size_type num_to_transport() const { return d_np_domain; }

    //! Total number of particles to transport in the entire problem/cycle.
    size_type total_num_to_transport() const { return d_np_total; }

    // >>> CLASS ACCESSORS

    //! Get the current fission site container.
    const Fission_Site_Container& fission_sites() const
    {
        if (!d_fission_sites)
            return d_dummy_container;
        return *d_fission_sites;
    }

    //! Total number of requested particles per cycle.
    size_type Np() const { return d_np_requested; }

    //! Number transported so far on this domain.
    size_type num_run() const { return d_num_run; }

    //! Number left to transport on this domain.
    size_type num_left() const { return d_num_left; }

    //! Get fission source lower coords for testing purposes
    const Space_Vector& lower_coords() const { return d_lower; }

    //! Get fission source width for testing purposes
    const Space_Vector& width() const { return d_width; }

    //! Set a new number per cycle.
    void update_Np(size_type np) { d_np_requested = np; }

  protected:
    // >>> IMPLEMENTATION

    typedef Source<Geometry> Base;
    using Base::b_geometry;
    using Base::b_physics;
    using Base::b_rng_control;
    using Base::b_nodes;
    using Base::make_RNG;

    // Build the domain replicated fission source.
    void build_DR(SP_Cart_Mesh mesh, Const_Array_View fis_dens);

    // Sample the geometry.
    int sample_geometry(Space_Vector &r, const Space_Vector &omega,
                        Particle_t &p, RNG_t rng);

    // Initial fission source lower coords and width.
    Space_Vector d_lower;
    Space_Vector d_width;

    // Requested particles per cycle.
    size_type d_np_requested;

    // Number of particles: total, domain
    size_type d_np_total, d_np_domain;

    // Particle weight.
    double d_wt;

    // Number of source particles left in the current domain.
    size_type d_num_left;

    // Number of particles run on the current domain.
    size_type d_num_run;

    // Dummy fission site container.
    Fission_Site_Container d_dummy_container;

    // Mesh-based starting distribution.
    int          d_current_cell;
    Vec_Int      d_fis_dist;
    SP_Cart_Mesh d_fis_mesh;
};

} // end namespace profugus

#endif // MC_mc_Fission_Source_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Source.hh
//---------------------------------------------------------------------------//
