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
#include <thrust/device_vector.h>

#include "Box_Shape.cuh"
#include "Source.cuh"
#include "Particle.cuh"
#include "Definitions.hh"

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
    typedef Source<Geometry>                     Base;
    typedef Geometry                             Geometry_t;
    typedef Teuchos::RCP<Teuchos::ParameterList> RCP_Std_DB;
    typedef Particle<Geometry>                   Particle_t;
    typedef cuda_utils::Space_Vector             Space_Vector;
    typedef std::shared_ptr<Box_Shape>           SP_Shape;
    typedef cuda::Shared_Device_Ptr<Box_Shape>   SDP_Shape;
    typedef cuda::Shared_Device_Ptr<Geometry_t>  SDP_Geometry;
    typedef std::shared_ptr<Particle_t>          SP_Particle;
    typedef def::Vec_Dbl                         Vec_Dbl;
    typedef def::size_type                       size_type;
    typedef thrust::device_vector<double>        Device_Dbl;
    //@}

  private:
    // >>> DATA

    // Source geometric shape.
    Box_Shape  *d_geo_shape;

    // Energy shape CDF.
    Device_Dbl d_erg_cdf_vec;
    double *d_erg_cdf;

  public:
    // Constructor.
    Uniform_Source(RCP_Std_DB db, SDP_Geometry geometry);

    // Build the initial source on the host.
    void build_source(SDP_Shape geometric_shape);

    // >>> REQUIRED PUBLIC INTERFACE

    // Get a particle from the source.
    __device__ inline Particle_t get_particle(
        std::size_t idx, RNG_State_t *rng) const;

    // >>> CLASS ACCESSORS

    //! Starting weight for histories
    __host__ __device__ double wt() const { return 1.0; }

  private:
    // >>> IMPLEMENTATION

    using Base::b_geometry;
    using Base::d_np_requested;
    using Base::d_np_total;
    using Base::d_np_domain;
    using Base::d_np_left;

    int d_num_groups;
};

} // end namespace cuda_mc

#include "Uniform_Source.i.cuh"

#endif // cuda_mc_Uniform_Source_cuh



//---------------------------------------------------------------------------//
//                 end of Uniform_Source.cuh
//---------------------------------------------------------------------------//
