//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Uniform_Source.cuh
 * \author Steven Hamilton
 * \date   Tue May 06 16:43:26 2014
 * \brief  Uniform_Source class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Uniform_Source_cuh
#define cuda_mc_Uniform_Source_cuh

#include <vector>

#include "Box_Shape.cuh"
#include "Source.cuh"
#include "Particle.hh"
#include "cuda_utils/Device_Vector.hh"

#include "Teuchos_ParameterList.hpp"

namespace cuda_mc
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
    typedef Source<Geometry>                    Base;
    typedef Geometry                            Geometry_t;
    typedef Teuchos::RCP<Teuchos::ParameterList> RCP_Std_DB;
    typedef Particle<Geometry>                  Particle_t;
    typedef typename Geometry_t::Space_Vector   Space_Vector;
    typedef std::shared_ptr<Box_Shape>          SP_Shape;
    typedef cuda::Shared_Device_Ptr<Box_Shape>  SDP_Shape;
    typedef cuda::Shared_Device_Ptr<Geometry_t> SDP_Geometry;
    typedef std::shared_ptr<Particle_t>         SP_Particle;
    typedef curandState_t                       RNG_t;
    typedef def::Vec_Dbl                        Vec_Dbl;
    typedef def::size_type                      size_type;
    typedef cuda::arch::Device                  Arch_t;
    typedef cuda::Device_Vector<Arch_t,double>  Device_Dbl;
    //@}

  private:
    // >>> DATA

    // Source geometric shape.
    SDP_Shape   d_geo_shape_host;
    Box_Shape  *d_geo_shape;

    // Energy shape CDF.
    Device_Dbl d_erg_cdf;

  public:
    // Constructor.
    Uniform_Source(RCP_Std_DB db, SDP_Geometry geometry);

    // Build the initial source on the host.
    void build_source(SP_Shape geometric_shape);

    // >>> REQUIRED PUBLIC INTERFACE

    // Get a particle from the source.
    __device__ Particle_t get_particle(RNG_t &rng);

    //! Boolean operator for source (true when source still has particles).
    __host__ __device__ bool empty() const { return d_np_left == 0; }

    //! Number of particles to transport in the source on the current domain.
    __host__ __device__
    size_type num_to_transport() const { return d_np_domain; }

    //! Total number of particles to transport in the entire problem/cycle.
    __host__ __device__
    size_type total_num_to_transport() const { return d_np_total; }

    // >>> CLASS ACCESSORS

    //! Total number of requested particles.
    __host__ __device__
    size_type Np() const { return d_np_requested; }

    //! Number transported so far on this domain.
    __host__ __device__
    size_type num_run() const { return d_np_run; }

    //! Number left to transport on this domain.
    __host__ __device__
    size_type num_left() const { return d_np_left; }

  private:
    // >>> IMPLEMENTATION
    int d_nodes;

    using Base::b_geometry;

    int d_num_groups;

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

} // end namespace cuda_mc

#include "Uniform_Source.i.cuh"

#endif // cuda_mc_Uniform_Source_cuh



//---------------------------------------------------------------------------//
//                 end of Uniform_Source.cuh
//---------------------------------------------------------------------------//
