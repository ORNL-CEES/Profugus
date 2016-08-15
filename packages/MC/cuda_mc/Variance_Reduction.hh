//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Variance_Reduction.hh
 * \author Thomas M. Evans
 * \date   Thu May 08 11:08:46 2014
 * \brief  Variance_Reduction base class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Variance_Reduction_hh
#define cuda_mc_Variance_Reduction_hh

#include <memory>

#include "utils/Definitions.hh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Shared_Device_Ptr.hh"
#include "cuda_utils/Stream.hh"
#include "Bank.hh"
#include "Physics.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Variance_Reduction
 * \brief Base class for applying variance reduction methods
 */
//===========================================================================//

template <class Geometry>
class Variance_Reduction
{
  public:
    //@{
    //! Useful typedefs.
    typedef Geometry                              Geometry_t;
    typedef typename Geometry_t::Geo_State_t      Geo_State_t;
    typedef Physics<Geometry>                     Physics_t;
    typedef typename Physics_t::Particle_Vector_t Particle_Vector_t;
    typedef Bank<Geometry>                        Bank_t;
    //@}

  public:

    // Virtual destructor.
    virtual ~Variance_Reduction() { /*...*/ }

    //! Set the geometry class
    void set( const cuda_utils::Shared_Device_Ptr<Geometry_t>& geometry) 
    { REQUIRE(geometry); b_geometry = geometry; }

    //! Set the physics class
    void set(const cuda_utils::Shared_Device_Ptr<Physics_t>& physics) 
    { REQUIRE(physics); b_physics = physics; }

    //! Return the splitting boolean
    bool uses_splitting() { return b_splitting; }

    // >>> VARIANCE REDUCTION INTERFACE

    //! Apply variance reduction method after a surface crossing
    virtual void post_surface(
	cuda_utils::Shared_Device_Ptr<Particle_Vector_t>& particles, 
	cuda_utils::Shared_Device_Ptr<Bank_t>& bank,
        cuda_utils::Stream<cuda_utils::arch::Device> stream =
        cuda_utils::Stream<cuda_utils::arch::Device>() ) const = 0;

    //! Apply variance reduction method after a collision
    virtual void post_collision(
	cuda_utils::Shared_Device_Ptr<Particle_Vector_t>& particles, 
	cuda_utils::Shared_Device_Ptr<Bank_t>& bank,
        cuda_utils::Stream<cuda_utils::arch::Device> stream =
        cuda_utils::Stream<cuda_utils::arch::Device>() ) const = 0;

  protected:
    // Problem geometry implementation.
    cuda_utils::Shared_Device_Ptr<Geometry_t> b_geometry;

    // Problem physics implementation.
    cuda_utils::Shared_Device_Ptr<Physics_t> b_physics;

    // Boolean is true if VR method uses splitting
    bool b_splitting;
};

} // end namespace cuda_profugus

#endif // cuda_mc_Variance_Reduction_hh

//---------------------------------------------------------------------------//
//                 end of Variance_Reduction.hh
//---------------------------------------------------------------------------//
