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
#include <thrust/device_vector.h>

#include "utils/View_Field.hh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Definitions.hh"
#include "cuda_utils/Device_Memory_Manager.hh"
#include "cuda_utils/Device_View_Field.hh"
#include "cuda_utils/Utility_Functions.hh"
#include "geometry/Definitions.hh"
#include "geometry/Bounding_Box.hh"
#include "Mesh_State.hh"
#include "Mesh_State_Vector.cuh"
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

    typedef size_t                                size_type;
    typedef profugus::geometry::cell_type         cell_type;
    typedef cuda_utils::Space_Vector              Space_Vector;
    typedef Mesh_State                            Geo_State_t;
    typedef Mesh_State_Vector                     Geo_State_Vector_t;
    typedef Mesh_State_Vector_DMM                 Geo_State_Vector_DMM_t;
    typedef cuda::const_Device_View_Field<int>    Int_View;
    typedef cuda::const_Device_View_Field<double> Double_View;

    //! Constructor
    Mesh_Geometry(Cartesian_Mesh mesh,
                  Int_View       matids,
                  Int_View       reflect)
      : d_mesh(mesh)
      , d_matids(matids)
      , d_reflect(reflect)
    {}

    // >>> DEVICE API

    //! Initialize track.
    __device__ inline void initialize(const Space_Vector&        r,
                                      const Space_Vector&        direction,
                                            Geo_State_Vector_t&  state_vector,
                                            int                  pid) const;

    //! Get distance to next boundary
    __device__ inline double distance_to_boundary(
        Geo_State_Vector_t& state_vector, int pid) const;

    //! Move to and cross a surface in the current direction.
    __device__ void move_to_surface(Geo_State_Vector_t& state_vector,
                                    int                 pid) const
    {
        move(state_vector.next_dist(pid), state_vector, pid);

        using def::I; using def::J; using def::K;

        const int num_cells_x = d_mesh.num_cells_along(I);
        const int num_cells_y = d_mesh.num_cells_along(J);
        const int num_cells_z = d_mesh.num_cells_along(K);

        constexpr int face_start = Geo_State_t::MINUS_X;
        auto& exiting_face    = state_vector.exiting_face(pid);
        auto& reflecting_face = state_vector.reflecting_face(pid);
        auto& next_ijk = state_vector.next_ijk(pid);
        exiting_face    = Geo_State_t::NONE;
        reflecting_face = Geo_State_t::NONE;
        if( next_ijk[I] < 0 )
        {
            exiting_face = Geo_State_t::MINUS_X;
            if( d_reflect[Geo_State_t::MINUS_X-face_start] )
                reflecting_face = Geo_State_t::MINUS_X;
        }
        else if( next_ijk[I] == num_cells_x )
        {
            exiting_face = Geo_State_t::PLUS_X;
            if( d_reflect[Geo_State_t::PLUS_X-face_start] )
                reflecting_face = Geo_State_t::PLUS_X;
        }
        else if( next_ijk[J] < 0 )
        {
            exiting_face = Geo_State_t::MINUS_Y;
            if( d_reflect[Geo_State_t::MINUS_Y-face_start] )
                reflecting_face = Geo_State_t::MINUS_Y;
        }
        else if( next_ijk[J] == num_cells_y )
        {
            exiting_face = Geo_State_t::PLUS_Y;
            if( d_reflect[Geo_State_t::PLUS_Y-face_start] )
                reflecting_face = Geo_State_t::PLUS_Y;
        }
        else if( next_ijk[K] < 0 )
        {
            exiting_face = Geo_State_t::MINUS_Z;
            if( d_reflect[Geo_State_t::MINUS_Z-face_start] )
                reflecting_face = Geo_State_t::MINUS_Z;
        }
        else if( next_ijk[K] == num_cells_z )
        {
            exiting_face = Geo_State_t::PLUS_Z;
            if( d_reflect[Geo_State_t::PLUS_Z-face_start] )
                reflecting_face = Geo_State_t::PLUS_Z;
        }

        // If we're not reflecting, update cell index
        if( reflecting_face == Geo_State_t::NONE )
            state_vector.ijk(pid) = state_vector.next_ijk(pid);
    }

    //! Move a distance \e d to a point in the current direction.
    __device__ void move_to_point(double              d,
                                  Geo_State_Vector_t& state_vector,
                                  int                 pid) const
    {
        move(d, state_vector, pid);

        update_state(state_vector, pid);
    }

    //! Number of cells (excluding "outside" cell)
    __device__ cell_type num_cells() const
    {
        return d_mesh.num_cells();
    }

    //! Return the current cell ID, valid only when inside the mesh
    __device__ cell_type cell(const Geo_State_Vector_t& state_vector,
                              int                       pid) const
    {
        DEVICE_REQUIRE(boundary_state(state_vector,pid) !=
                       profugus::geometry::OUTSIDE);

        using def::I; using def::J; using def::K;
        cell_type c = num_cells();
        const auto& ijk = state_vector.ijk(pid);
        bool found = d_mesh.index(ijk[I], ijk[J], ijk[K], c);

        DEVICE_ENSURE(found);
        return c;
    }

    //! Return the current material ID
    __device__ int matid(const Geo_State_Vector_t& state_vector, int pid) const
    {
        DEVICE_REQUIRE(cell(state_vector,pid) < num_cells());

        return d_matids[cell(state_vector,pid)];
    }

    //! Return the state with respect to outer geometry boundary
    __device__ profugus::geometry::Boundary_State boundary_state(
        const Geo_State_Vector_t& state_vector, int pid) const
    {
        using def::I; using def::J; using def::K;

        const auto& ijk = state_vector.ijk(pid);
        if (state_vector.reflecting_face(pid) != Geo_State_t::NONE)
        {
            return profugus::geometry::REFLECT;
        }
        else if (   (ijk[I] == -1)
                 || (ijk[J] == -1)
                 || (ijk[K] == -1)
                 || (ijk[I] == d_mesh.num_cells_along(I))
                 || (ijk[J] == d_mesh.num_cells_along(J))
                 || (ijk[K] == d_mesh.num_cells_along(K)))
        {
            return profugus::geometry::OUTSIDE;
        }
        return profugus::geometry::INSIDE;
    }

    //! Return the current position.
    __device__ Space_Vector position(const Geo_State_Vector_t& state_vector,
                                     int                       pid) const
    {
        return state_vector.pos(pid);
    }

    //! Return the current direction.
    __device__ Space_Vector direction(const Geo_State_Vector_t& state_vector,
                                      int                       pid) const
    {
        return state_vector.dir(pid);
    }

    //! Change the direction to \p new_direction.
    __device__ void change_direction( const Space_Vector&        new_direction,
                                            Geo_State_Vector_t&  state_vector,
                                            int                  pid) const
    {
        // update and normalize the direction
        state_vector.dir(pid) = new_direction;
        cuda::utility::vector_normalize(state_vector.dir(pid));
    }

    //! Change the direction through an angle
    __device__ void change_direction( double               costheta,
                                      double               phi,
                                      Geo_State_Vector_t&  state_vector,
                                      int                  pid) const
    {
        cuda::utility::cartesian_vector_transform(costheta, phi,
                                                  state_vector.dir(pid));
    }

    //! Reflect the direction at a reflecting surface.
    __device__ inline bool reflect(Geo_State_Vector_t& state_vector,
                                   int                 pid) const;

    //! Return the outward normal at the location dictated by the state.
    __device__ inline Space_Vector normal(const Geo_State_Vector_t& state_vector,
                                          int                       pid) const;

  private:

    // >>> IMPLEMENTATION

    // Update state tracking information
    __device__ void update_state(Geo_State_Vector_t &state_vector,
                                 int                 pid) const
    {
        // Find the logical indices of the cell along each axis
        d_mesh.find_upper(state_vector.pos(pid), state_vector.ijk(pid));
    }

    //! Move a particle a distance \e d in the current direction.
    __device__ void move(double              dist,
                         Geo_State_Vector_t &state_vector,
                         int                 pid) const
    {
        using def::I; using def::J; using def::K;

        auto &pos       = state_vector.pos(pid);
        const auto &dir = state_vector.dir(pid);

        DEVICE_REQUIRE(dist >= 0.0);
        DEVICE_REQUIRE(cuda::utility::soft_equiv(
                       cuda::utility::vector_magnitude(dir), 1.0, 1.0e-6));

        // advance the particle
        for (int dim : {I, J, K})
            pos[dim] += dist * dir[dim];
    }

    Cartesian_Mesh d_mesh;

    Int_View d_matids;
    Int_View d_reflect;
};

//===========================================================================//
/*!
 * \class Mesh_Geometry_DMM
 * \brief Device memory manager for Mesh_Geometry
 */
//===========================================================================//

class Mesh_Geometry_DMM : public cuda::Device_Memory_Manager<Mesh_Geometry>
{
  public:

    typedef Mesh_Geometry                        Geometry_t;
    typedef size_t                               size_type;
    typedef profugus::geometry::cell_type        cell_type;
    typedef std::vector<double>                  Vec_Dbl;
    typedef std::vector<int>                     Vec_Int;
    typedef cuda::arch::Device                   Arch;
    typedef thrust::device_vector<double>        Dev_Dbl_Vec;
    typedef thrust::device_vector<int>           Dev_Int_Vec;

    //! Constructor
    Mesh_Geometry_DMM(const Vec_Dbl &x_edges,
                      const Vec_Dbl &y_edges,
                      const Vec_Dbl &z_edges);

    // DMM Interface
    Mesh_Geometry device_instance()
    {
        Mesh_Geometry geom(d_mesh.device_instance(),
                           cuda::make_view(d_matids),
                           cuda::make_view(d_reflect));
        return geom;
    }

    // Number of mesh cells in geometry
    size_type num_cells() const
    {
        return d_mesh.num_cells();
    }

    //! Set materials
    void set_matids(const Vec_Int& matids)
    {
        DEVICE_REQUIRE(matids.size() == d_mesh.num_cells());
        d_matids = matids;
    }

    //! Set reflecting boundaries
    void set_reflecting(const Vec_Int &reflecting_faces)
    {
        DEVICE_REQUIRE( reflecting_faces.size() == 6 );
        d_reflect = reflecting_faces;
    }

    //! Access the underlying mesh directly
    const Cartesian_Mesh_DMM& mesh() const {return d_mesh;}

    //! Low corner of problem domain
    def::Space_Vector lower() const
    {
        return d_mesh.lower();
    }

    //! High corner of problem domain
    def::Space_Vector upper() const
    {
        return d_mesh.upper();
    }

    //! All cell volumes
    const std::vector<double>& volumes() const {return d_mesh.volumes();}

    Cartesian_Mesh_DMM d_mesh;

    thrust::device_vector<int> d_matids;
    thrust::device_vector<int> d_reflect;
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
