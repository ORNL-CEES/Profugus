//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Uniform_Source.hh
 * \author Thomas M. Evans
 * \date   Tue May 06 16:43:26 2014
 * \brief  Uniform_Source class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Uniform_Source_hh
#define mc_Uniform_Source_hh

#include <vector>

#include "Shape.hh"
#include "Source.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Uniform_Source
 * \brief Defines a uniformly sampled source for fixed-source problems.
 *
 * Currently, the only implemented source shape is a box (see Box_Shape).
 * Also, all angular sampling is isotropic.
 *
 * \section uniform_source_db Standard DB Entries for Uniform_Source
 *
 * The following standard data entries control the Uniform_Source:
 *
 * \arg \c Np (int) number of particles to use in each cycle (default:
 * 1000)
 *
 * \arg \c spectral_shape (Array<double>) source spectral (energy) shape by
 * group (default: flat)
 */
/*!
 * \example mc/test/tstUniform_Source.cc
 *
 * Test of Uniform_Source.
 */
//===========================================================================//

template <class Geometry>
class Uniform_Source : public Source<Geometry>
{
  public:
    //@{
    //! Typedefs.
    typedef Geometry                            Geometry_t;
    typedef Physics<Geometry_t>                 Physics_t;
    typedef typename Physics_t::RCP_Std_DB      RCP_Std_DB;
    typedef typename Physics_t::Particle_t      Particle_t;
    typedef typename Geometry_t::Space_Vector   Space_Vector;
    typedef std::shared_ptr<Shape>              SP_Shape;
    typedef std::shared_ptr<Geometry_t>         SP_Geometry;
    typedef std::shared_ptr<Physics_t>          SP_Physics;
    typedef std::shared_ptr<RNG_Control>        SP_RNG_Control;
    typedef std::shared_ptr<Particle_t>         SP_Particle;
    typedef def::Vec_Dbl                        Vec_Dbl;
    typedef def::size_type                      size_type;
    //@}

  private:
    // >>> DATA

    // Source geometric shape.
    SP_Shape d_geo_shape;

    // Energy shape CDF.
    Vec_Dbl d_erg_cdf;

  public:
    // Constructor.
    Uniform_Source(RCP_Std_DB db, SP_Geometry geometry, SP_Physics physics,
                   SP_RNG_Control rng_control);

    // Build the initial source.
    void build_source(SP_Shape geometric_shape);

    // >>> DERIVED PUBLIC INTERFACE

    // Get a particle from the source.
    SP_Particle get_particle();

    //! Boolean operator for source (true when source still has particles).
    bool empty() const { return d_np_left == 0; }

    //! Number of particles to transport in the source on the current domain.
    size_type num_to_transport() const { return d_np_domain; }

    //! Total number of particles to transport in the entire problem/cycle.
    size_type total_num_to_transport() const { return d_np_total; }

    // >>> CLASS ACCESSORS

    //! Total number of requested particles.
    size_type Np() const { return d_np_requested; }

    //! Number transported so far on this domain.
    size_type num_run() const { return d_np_run; }

    //! Number left to transport on this domain.
    size_type num_left() const { return d_np_left; }

  private:
    // >>> IMPLEMENTATION

    typedef Source<Geometry> Base;
    using Base::b_geometry;
    using Base::b_physics;
    using Base::b_nodes;

    // Build the domain replicated source.
    void build_DR();

    // Requested particles.
    size_type d_np_requested;

    // Number of particles: total, domain
    size_type d_np_total;
    size_type d_np_domain;

    // Particle weight.
    double d_wt;

    // Number of source particles left in the current domain.
    size_type d_np_left;

    // Number of particles run on the current domain.
    size_type d_np_run;
};

} // end namespace profugus

#endif // mc_Uniform_Source_hh

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.hh
//---------------------------------------------------------------------------//
