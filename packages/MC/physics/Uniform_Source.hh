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
#include <atomic>

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

class Uniform_Source : public Source
{
  public:
    //@{
    //! Typedefs.
    typedef Physics_t::RCP_Std_DB  RCP_Std_DB;
    typedef def::Vec_Dbl           Vec_Dbl;
    typedef std::shared_ptr<Shape> SP_Shape;
    //@}

  private:
    // >>> DATA

    // Source geometric shape.
    SP_Shape d_geo_shape;

    // Energy shape CDF.
    Vec_Dbl d_erg_cdf;

  public:
    // Constructor.
    Uniform_Source(RCP_Std_DB db, SP_Geometry geometry, SP_Physics physics);

    // Build the initial source.
    void build_source(SP_Shape geometric_shape);

    // >>> DERIVED PUBLIC INTERFACE

    // Get a particle from the source.
    SP_Particle get_particle();

    //! Boolean operator for source (true when source still has particles).
    bool empty() const { return d_np_run.load() == d_np_domain; }

    //! Number of particles to transport in the source on the current domain.
    size_type num_to_transport() const { return d_np_domain; }

    //! Total number of particles to transport in the entire problem/cycle.
    size_type total_num_to_transport() const { return d_np_total; }

    // >>> CLASS ACCESSORS

    //! Total number of requested particles.
    size_type Np() const { return d_np_requested; }

    //! Number transported so far on this domain.
    size_type num_run() const { return d_np_run.load(); }

  private:
    // >>> IMPLEMENTATION

    typedef Source Base;

    // Build the domain replicated source.
    void build_DR();

    // Requested particles.
    size_type d_np_requested;

    // Number of particles: total, domain
    size_type d_np_total;
    size_type d_np_domain;

    // Particle weight.
    const double d_wt;

    // Number of particles run on the current domain.
    std::atomic<size_type> d_np_run;
};

} // end namespace profugus

#endif // mc_Uniform_Source_hh

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.hh
//---------------------------------------------------------------------------//
