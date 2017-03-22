//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_Geometry.cuh
 * \author Tom Evans
 * \date   Mon Jan 30 00:14:27 2017
 * \brief  RTK_Geometry CUDA-device class declarations.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_RTK_Geometry_cuh
#define MC_cuda_rtk_RTK_Geometry_cuh

#include "MC/geometry/RTK_Geometry.hh"
#include "RTK_Array.cuh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class RTK_Geometry
 * \brief CUDA implementation of profugus::RTK_Geometry.
 *
 * For cuda implementations the RTK_Geometry is explicity built using
 * Core_Array.  Thus, we \b always have a 2-level geometry, which is
 * sufficient for PWR core-level problems and BWR lattices.
 */
/*!
 * \example cuda_rtk/test/tstRTK_Geometry_cuda.cc
 *
 * Test of RTK_Geometry.
 */
//===========================================================================//

class RTK_Geometry
{
  public:
    //@{
    using Array_t        = Core_Array;
    using Geo_State_t    = RTK_State;
    using Space_Vector   = Geo_State_t::Space_Vector;
    using Boundary_State = profugus::geometry::Boundary_State;
    //@}

  private:
    // >>> DATA

    // Underlying core array.
    Array_t d_array;

    // Level of array.
    const int d_level;

  public:
    // Constructor.
    explicit RTK_Geometry(Array_t array);

    // >>> DEVICE FUNCTIONS

    // GEOMETRY INTERFACE FUNCTIONS

    // Initialize a track.
    __device__
    inline void initialize(const Space_Vector &r, const Space_Vector &direction,
                           Geo_State_t &state) const;

    //! Get distance to next boundary.
    __device__
    double distance_to_boundary(Geo_State_t &state) const
    {
        d_array.distance_to_boundary(state.d_r, state.d_dir, state);
        return state.dist_to_next_region;
    }

    //! Move to and cross a cell surface.
    __device__
    void move_to_surface(Geo_State_t &state) const
    {
        // move the particle
        move(state.dist_to_next_region, state);

        // process the particle through the surface
        d_array.cross_surface(state.d_r, state);
    }

    //! Move the particle to a point in the current direction.
    __device__
    void move_to_point(double d, Geo_State_t &state) const
    {
        // move the particle
        move(d, state);

        // update the array state to clear any surface tags
        d_array.update_state(state);
    }

    //! Return the current material ID.
    __device__
    int matid(const Geo_State_t &state) const
    {
        return d_array.matid(state);
    }

    //! Return the current cell ID.
    __device__
    int cell(const Geo_State_t &state) const
    {
        return d_array.cellid(state);
    }

    // Return the state with respect to outer geometry boundary.
    __device__
    inline Boundary_State boundary_state(const Geo_State_t &state) const;

    //! Return the current position.
    __device__
    Space_Vector position(const Geo_State_t &state) const { return state.d_r; }

    //! Return the current direction.
    __device__
    Space_Vector direction(const Geo_State_t &state) const {return state.d_dir;}

    //! Change the particle direction.
    __device__
    void change_direction(const Space_Vector &new_direction,
                                Geo_State_t  &state) const
    {
        // update the direction
        state.d_dir = new_direction;

        // normalize the direction
        cuda_utils::utility::vector_normalize(state.d_dir);
    }

    // Change the direction through angles \f$(\theta,\phi)\f$.
    __device__
    void change_direction(double costheta, double phi, Geo_State_t &state) const
    {
        cuda_utils::utility::cartesian_vector_transform(costheta, phi, state.d_dir);
    }

    // Reflect the direction at a reflecting surface.
    __device__
    inline bool reflect(Geo_State_t &state) const;

    // Return the outward normal.
    __device__
    inline Space_Vector normal(const Geo_State_t &state) const;

  private:
    // >>> IMPLEMENTATION

    //! Move a particle a distance \e d in the current direction.
    __device__
    void move(double d, Geo_State_t &state) const
    {
        DEVICE_REQUIRE(d >= 0.0);
        DEVICE_REQUIRE(cuda_utils::utility::soft_equiv(
                       cuda_utils::utility::vector_magnitude(state.d_dir),
                       1.0, 1.0e-6));

        // advance the particle (unrolled loop)
        state.d_r[def::X] += d * state.d_dir[def::X];
        state.d_r[def::Y] += d * state.d_dir[def::Y];
        state.d_r[def::Z] += d * state.d_dir[def::Z];
    }
};

//===========================================================================//
/*!
 * \class RTK_Geometry_DMM
 * \brief Device_Memory_Manager for an RTK_Geometry
 */
//===========================================================================//

class RTK_Geometry_DMM
    : public cuda::Device_Memory_Manager<RTK_Geometry>
{
  public:
    using Geometry_t    = RTK_Geometry;

  private:
    // Types.
    using Base          = cuda::Device_Memory_Manager<RTK_Geometry>;
    using Array_DMM     = Core_Array_DMM;
    using Host_Geometry = profugus::Core;

  private:
    // >>> DEVICE MEMORY MANAGEMENT DATA

    // Memory management of underlying array
    Array_DMM d_array_manager;

    // Vector of all cell volumes
    std::vector<double> d_volumes;

    // Extents
    def::Space_Vector d_lower;
    def::Space_Vector d_upper;

  public:
    // Constructor.
    explicit RTK_Geometry_DMM(const Host_Geometry &geometry);

    // >>> DERIVED INTERFACE

    // Create a device instance.
    RTK_Geometry device_instance();

    const std::vector<double>& volumes() const
    {
        return d_volumes;
    }

    // Geometry extents
    def::Space_Vector lower() const {return d_lower;}
    def::Space_Vector upper() const {return d_upper;}
};

//---------------------------------------------------------------------------//
// SUPPORTED GEOMETRIES
//---------------------------------------------------------------------------//

//! Geometry built on Core_Array.
using Core = RTK_Geometry;

//! Geometry memory manager for Core.
using Core_DMM = RTK_Geometry_DMM;

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// CUDA FUNCTIONS
//---------------------------------------------------------------------------//

#include "RTK_Geometry.i.cuh"

#endif // MC_cuda_rtk_RTK_Geometry_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Geometry.cuh
//---------------------------------------------------------------------------//
