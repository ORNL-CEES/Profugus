//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/Mesh_Geometry.hh
 * \author Steven Hamilton
 * \date   Tue Dec 15 14:10:10 2015
 * \brief  Mesh_Geometry class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_geometry_Mesh_Geometry_hh
#define MC_cuda_geometry_Mesh_Geometry_hh

#include <vector>

#include "utils/View_Field.hh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Definitions.hh"
#include "cuda_utils/Device_Vector.hh"
#include "cuda_utils/Utility_Functions.hh"
#include "geometry/Definitions.hh"
#include "geometry/Bounding_Box.hh"
#include "Mesh_State.hh"
#include "Cartesian_Mesh.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Mesh_Geometry
 * \brief Track particles through structured Cartesian Mesh.
 *
 * \sa Mesh_Geometry.cu for detailed descriptions.
 */
/*!
 * \example cuda_geometry/test/tstMesh_Geometry_cuda.cc
 *
 * Test of Mesh_Geometry.
 */
//===========================================================================//

class Mesh_Geometry
{
  public:

    typedef size_t                               size_type;
    typedef profugus::geometry::cell_type        cell_type;
    typedef std::vector<double>                  Vec_Dbl;
    typedef std::vector<int>                     Vec_Int;
    typedef cuda::arch::Device                   Arch;
    typedef cuda::Device_Vector<Arch,double>     Dev_Dbl_Vec;
    typedef std::shared_ptr<Dev_Dbl_Vec>         SP_Dev_Dbl_Vec;
    typedef cuda::Device_Vector<Arch,int>        Dev_Int_Vec;
    typedef std::shared_ptr<Dev_Int_Vec>         SP_Dev_Int_Vec;
    typedef cuda::Coordinates                    Coordinates;
    typedef cuda::Space_Vector                   Space_Vector;
    typedef Mesh_State                           Geo_State_t;


    //! Constructor
    Mesh_Geometry(const Vec_Dbl &x_edges,
                  const Vec_Dbl &y_edges,
                  const Vec_Dbl &z_edges);

    Mesh_Geometry(const Mesh_Geometry &geom)             = default;
    Mesh_Geometry& operator=( const Mesh_Geometry &geom) = default;

    // >>> HOST API

    //! Set materials
    void set_matids(const Vec_Int & matids)
    {
        REQUIRE( matids.size() == num_cells() );
        d_matid_vec =
            std::make_shared<Dev_Int_Vec>(profugus::make_view(matids));
        dd_matids = d_matid_vec->data();
    }

    //! Access the underlying mesh directly
    const Cartesian_Mesh& mesh() const { return d_mesh; }

    //! Low corner of problem domain
    cuda::Space_Vector lower() const
    {
        return d_mesh.lower();
    }

    //! High corner of problem domain
    cuda::Space_Vector upper() const
    {
        return d_mesh.upper();
    }

    // >>> DEVICE API

    //! Initialize track.
    __device__ inline void initialize(const Space_Vector& r,
                                      const Space_Vector& direction,
                                            Geo_State_t&  state) const;

    //! Get distance to next boundary
    __device__ inline double distance_to_boundary(Geo_State_t& state) const;

    //! Move to and cross a surface in the current direction.
    __device__ void move_to_surface(Geo_State_t& state) const
    {
        move(state.next_dist, state);
        state.ijk = state.next_ijk;
    }

    //! Move a distance \e d to a point in the current direction.
    __device__ void move_to_point(double d, Geo_State_t& state) const
    {
        move(d, state);

        update_state(state);
    }

    //! Number of cells (excluding "outside" cell)
    __host__ __device__ cell_type num_cells() const
    {
        return d_mesh.num_cells();
    }

    //! Return the current cell ID, valid only when inside the mesh
    __device__ cell_type cell(const Geo_State_t& state) const
    {
        REQUIRE(boundary_state(state) != profugus::geometry::OUTSIDE);

        using def::I; using def::J; using def::K;
        cell_type c = num_cells();
        bool found = d_mesh.index(state.ijk.i, state.ijk.j, state.ijk.k, c);

        ENSURE(found);
        return c;
    }

    //! Return the current material ID
    __device__ int matid(const Geo_State_t& state) const
    {
        REQUIRE(cell(state) < num_cells());

        return dd_matids[cell(state)];
    }

    //! Return the state with respect to outer geometry boundary
    __device__ profugus::geometry::Boundary_State boundary_state(
        const Geo_State_t& state) const
    {
        using def::I; using def::J; using def::K;

        if (   (state.ijk.i == -1)
            || (state.ijk.j == -1)
            || (state.ijk.k == -1)
            || (state.ijk.i == d_mesh.num_cells_along(I))
            || (state.ijk.j == d_mesh.num_cells_along(J))
            || (state.ijk.k == d_mesh.num_cells_along(K)))
        {
            return profugus::geometry::OUTSIDE;
        }
        return profugus::geometry::INSIDE;
    }

    //! Return the current position.
    __device__ Space_Vector position(const Geo_State_t& state) const
    {
        return state.d_r;
    }

    //! Return the current direction.
    __device__ Space_Vector direction(const Geo_State_t& state) const
    {
        return state.d_dir;
    }

    //! Change the direction to \p new_direction.
    __device__ void change_direction( const Space_Vector& new_direction,
                                            Geo_State_t&  state) const
    {
        // update and normalize the direction
        state.d_dir = new_direction;
        cuda::utility::vector_normalize(state.d_dir);
    }

    //! Change the direction through an angle
    __device__ void change_direction( double       costheta,
                                      double       phi,
                                      Geo_State_t& state) const
    {
        cuda::utility::cartesian_vector_transform(costheta, phi, state.d_dir);
    }

    //! Reflect the direction at a reflecting surface.
    __device__ bool reflect(Geo_State_t& state) const
    {
        INSIST(false, "no reflect in mesh_geometry");
        return false;
    }

    //! Return the outward normal at the location dictated by the state.
    __device__ Space_Vector normal(const Geo_State_t& state) const
    {
        INSIST(false, "normal in mesh_geometry");
        Space_Vector dir = {0.0, 0.0, 0.0};
        return dir;
    }

    // If the particle is outside the geometry, find distance
    __device__ double distance_to_interior(Geo_State_t &state);

    //! Get cell volume
    __device__ double volume(size_type cell) const
    {
        return d_mesh.volume(cell);
    }

  private:

    // >>> IMPLEMENTATION

    // Update state tracking information
    __device__ void update_state(Geo_State_t &state) const
    {
        using def::I; using def::J; using def::K;

        // Find the logical indices of the cell along each axis
        d_mesh.find_upper(state.d_r, state.ijk);
    }

    //! Move a particle a distance \e d in the current direction.
    __device__ void move(double dist, Geo_State_t &state) const
    {
        REQUIRE(dist >= 0.0);
        REQUIRE(cuda::utility::soft_equiv(
                cuda::utility::vector_magnitude(state.d_dir),
                1.0, 1.0e-6));

        // advance the particle (unrolled loop)
        state.d_r.x += dist * state.d_dir.x;
        state.d_r.y += dist * state.d_dir.y;
        state.d_r.z += dist * state.d_dir.z;
    }

    Cartesian_Mesh d_mesh;

    SP_Dev_Int_Vec d_matid_vec;
    int *dd_matids;
};

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
#include "Mesh_Geometry.i.hh"
//---------------------------------------------------------------------------//
#endif // MC_cuda_geometry_Mesh_Geometry_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/Mesh_Geometry.hh
//---------------------------------------------------------------------------//
