//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Uniform_Source.hh
 * \author Stuart Slattery
 * \brief  Uniform_Source class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Uniform_Source_hh
#define cuda_mc_Uniform_Source_hh

#include "utils/Definitions.hh"

#include "cuda_utils/Definitions.hh"
#include "cuda_utils/Shared_Device_Ptr.hh"
#include "cuda_utils/Stream.hh"

#include "Source.hh"
#include "Definitions.hh"
#include "Particle_Vector.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace cuda_profugus
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
template <class Geometry, class Shape>
class Uniform_Source : public Source<Geometry>
{
  public:
    //@{
    //! Typedefs.
    typedef Geometry                            Geometry_t;
    typedef Teuchos::ParameterList              ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t>       RCP_Std_DB;
    typedef typename Geometry_t::Space_Vector   Space_Vector;
    typedef events::Event                       Event_t;
    //@}

  private:
    // >>> DATA

    // Geometry.
    cuda_utils::Shared_Device_Ptr<Geometry> d_geometry;

    // Source geometric shape.
    cuda_utils::Shared_Device_Ptr<Shape> d_shape;

    // Number of energy groups.
    int d_num_groups;

    // Number of particle batches.
    int d_num_batch;

    // Particle batch size.
    int d_batch_size;

    // Energy shape CDF. On-device.
    double* d_erg_cdf;

  public:
    
    // Constructor.
    Uniform_Source( const RCP_Std_DB& db, 
		    const cuda_utils::Shared_Device_Ptr<Geometry>& geometry,
		    const int num_groups,
		    const int num_batch );

    // Destructor.
    ~Uniform_Source();

    // Build the initial source.
    void build_source( const cuda_utils::Shared_Device_Ptr<Shape>& shape );

    // >>> HOST-API

    // Get particles from the source.
    void get_particles( 
	cuda_utils::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles ) override;

    //! Boolean operator for source (true when source still has particles).
    bool empty() const override { return d_np_left == 0; }

    //! Number of particle batches.
    int num_batch() const { return d_num_batch; }

    //! Number of particles to transport in the source on the current domain.
    int num_to_transport() const override { return d_np_domain; }

    //! Total number of particles to transport in the entire problem/cycle.
    int total_num_to_transport() const override { return d_np_total; }

    // >>> CLASS ACCESSORS

    //! Total number of requested particles.
    int Np() const override { return d_np_requested; }

    //! Number transported so far on this domain.
    int num_run() const override { return d_np_run; }

    //! Number left to transport on this domain.
    int num_left() const override { return d_np_left; }

  private:
    // >>> IMPLEMENTATION

    // Build the domain replicated source.
    void build_DR();

  private:

    // Execution stream.
    cuda_utils::Stream<cuda_utils::arch::Device> d_stream;

    // Requested particles.
    int d_np_requested;

    // Number of particles: total, domain
    int d_np_total;
    int d_np_domain;

    // Particle weight.
    double d_wt;

    // Number of source particles left in the current domain.
    int d_np_left;

    // Number of particles run on the current domain.
    int d_np_run;
};

} // end namespace cuda_profugus

#endif // cuda_mc_Uniform_Source_hh

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.hh
//---------------------------------------------------------------------------//
