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

#include "cuda_utils/Device_Memory_Manager.hh"
#include "cuda_utils/Device_View_Field.hh"
#include "Box_Shape.cuh"
#include "Source.cuh"
#include "Particle_Vector.cuh"
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
class Uniform_Source
{
  public:
    //@{
    //! Typedefs.
    typedef Geometry                              Geometry_t;
    typedef Particle_Vector<Geometry>             Particle_Vector_t;
    typedef cuda_utils::Space_Vector              Space_Vector;
    typedef def::size_type                        size_type;
    typedef cuda::const_Device_View_Field<double> Double_View;
    //@}

  private:
    // >>> DATA

    // Geometry
    const Geometry_t *d_geometry;

    // Source geometric shape.
    Box_Shape  *d_geo_shape;

    // Energy shape CDF.
    Double_View d_erg_cdf;

  public:
    // Constructor.
    Uniform_Source(const Geometry_t *geometry,
                   Box_Shape        *geo_shape,
                   Double_View       erg_cdf)
      : d_geometry(geometry)
      , d_geo_shape(geo_shape)
      , d_erg_cdf(erg_cdf)
    {}

    // >>> REQUIRED PUBLIC INTERFACE

    // Get a particle from the source.
    __device__ inline void build_particle(
        int idx, Particle_Vector_t *particles) const;

    // >>> CLASS ACCESSORS

    //! Starting weight for histories
    __device__ double wt() const { return 1.0; }
};

//===========================================================================//
/*!
 * \class Uniform_Source_DMM
 * \brief Device memory manager for Uniform_Source.
 */
//===========================================================================//

template <class Geometry>
class Uniform_Source_DMM : public Source<Geometry>,
    public cuda::Device_Memory_Manager<Uniform_Source<Geometry>>
{
  public:
    //@{
    //! Typedefs.
    typedef Source<Geometry>                     Base;
    typedef Geometry                             Geometry_t;
    typedef Teuchos::RCP<Teuchos::ParameterList> RCP_Std_DB;
    typedef cuda_utils::Space_Vector             Space_Vector;
    typedef std::shared_ptr<Box_Shape>           SP_Shape;
    typedef cuda::Shared_Device_Ptr<Box_Shape>   SDP_Shape;
    typedef cuda::Shared_Device_Ptr<Geometry_t>  SDP_Geometry;
    typedef def::Vec_Dbl                         Vec_Dbl;
    typedef def::size_type                       size_type;
    typedef thrust::device_vector<double>        Device_Dbl;
    //@}

  private:
    // >>> DATA

    // Source geometric shape.
    SDP_Shape d_geo_shape;

    // Energy shape CDF.
    Device_Dbl d_erg_cdf;

  public:
    // Constructor.
    Uniform_Source_DMM(RCP_Std_DB db, SDP_Geometry geometry);

    // DMM interface
    Uniform_Source<Geometry_t> device_instance()
    {
        Uniform_Source<Geometry_t> src(b_geometry.get_device_ptr(),
                                       d_geo_shape.get_device_ptr(),
                                       cuda::make_view(d_erg_cdf));
        return src;
    }

    // Build the initial source on the host.
    void build_source(SDP_Shape geometric_shape);

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
