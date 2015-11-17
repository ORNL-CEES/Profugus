//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/Mesh_Geometry.hh
 * \author Thomas M. Evans
 * \date   Mon Jul 21 17:56:40 2014
 * \brief  Mesh_Geometry class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_Mesh_Geometry_hh
#define geometry_Mesh_Geometry_hh

#include <vector>
#include <memory>

#include "harness/DBC.hh"
#include "utils/Vector_Lite.hh"
#include "utils/Vector_Functions.hh"
#include "Cartesian_Mesh.hh"
#include "Mesh_State.hh"
#include "Tracking_Geometry.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Mesh_Geometry
 * \brief Track particles through a structured Cartesian mesh
 *
 * This class is designed to be used for mesh tally tracking,
 * material/cell discretization, etc.
 */
/*!
 * \example geometry/test/tstMesh_Geometry.cc
 *
 * Test of Mesh_Geometry.
 */
//===========================================================================//

class Mesh_Geometry : public Tracking_Geometry<Mesh_State>
{
  public:
    //@{
    //! Typedefs
    typedef def::Vec_Dbl             Vec_Dbl;
    typedef def::Vec_Int             Vec_Int;
    typedef std::shared_ptr<Vec_Int> SP_Vec_Int;
    typedef std::shared_ptr<Vec_Dbl> SP_Vec_Dbl;
    //@}

  private:
    // >>> DATA

    // Underlying mesh
    Cartesian_Mesh d_mesh;

    // Material specification (optional)
    SP_Vec_Int d_materials;

    // Cell volumes.
    SP_Vec_Dbl d_volumes;

  public:
    //! Construct with global edges
    Mesh_Geometry(const Vec_Dbl& x_edges,
                  const Vec_Dbl& y_edges,
                  const Vec_Dbl& z_edges);

    // Set materials (optional).
    void set_matids(SP_Vec_Int matids);

    // >>> DERIVED INTERFACE from Tracking_Geometry

    //! Initialize track.
    void initialize(const Space_Vector& r,
                    const Space_Vector& direction,
                    Geo_State_t       & state) const;

    //! Get distance to next boundary.
    double distance_to_boundary(Geo_State_t& state);

    //! Move to and cross a surface in the current direction.
    void move_to_surface(Geo_State_t& state)
    {
        move(state.next_dist, state);
        state.ijk = state.next_ijk;
    }

    //! Move a distance \e d to a point in the current direction.
    void move_to_point(double d, Geo_State_t& state)
    {
        move(d, state);

        update_state(state);
    }

    //! Number of cells (excluding "outside" cell)
    geometry::cell_type num_cells() const { return d_mesh.num_cells(); }

    //! Return the current cell ID, valid only when inside the mesh
    geometry::cell_type cell(const Geo_State_t &state) const
    {
        REQUIRE(boundary_state(state) != geometry::OUTSIDE);

        using def::I; using def::J; using def::K;
        Cartesian_Mesh::size_type c = num_cells();
        bool found = d_mesh.index(state.ijk[I], state.ijk[J], state.ijk[K], c);

        ENSURE(found);
        return c;
    }

    //! Return the current cell ID from a position inside the mesh
    geometry::cell_type cell(const Space_Vector &r) const
    {
        return Tracking_Geometry<Mesh_State>::cell(r);
    }

    //! Return the current material ID
    geometry::matid_type matid(const Geo_State_t &state) const
    {
        INSIST(d_materials, "Material IDs haven't been assigned");
        REQUIRE(cell(state) < d_materials->size());

        ENSURE((*d_materials)[cell(state)] >= 0);
        return (*d_materials)[cell(state)];
    }

    //! Return the material ID for the given location
    geometry::matid_type matid(const Space_Vector &r) const
    {
        return Tracking_Geometry<Mesh_State>::matid(r);
    }

    //! Return the state with respect to outer geometry boundary
    geometry::Boundary_State boundary_state(const Geo_State_t &state) const
    {
        using def::I; using def::J; using def::K;
        const Cartesian_Mesh::Dim_Vector& extents = d_mesh.extents();

        if (   (state.ijk[I] == -1)
            || (state.ijk[J] == -1)
            || (state.ijk[K] == -1)
            || (state.ijk[I] == extents[I])
            || (state.ijk[J] == extents[J])
            || (state.ijk[K] == extents[K]))
        {
            return geometry::OUTSIDE;
        }
        return geometry::INSIDE;
    }

    //! Return the boundary state for the given location
    geometry::Boundary_State boundary_state(const Space_Vector &r) const
    {
        return Tracking_Geometry<Mesh_State>::boundary_state(r);
    }

    //! Return the current position.
    Space_Vector position(const Geo_State_t& state) const
    {
        return state.pos;
    }

    //! Return the current direction.
    Space_Vector direction(const Geo_State_t& state) const
    {
        return state.dir;
    }

    //! Change the direction to \p new_direction.
    void change_direction(
            const Space_Vector& new_direction,
            Geo_State_t& state)
    {
        // update and mnormalizethe direction
        state.dir = new_direction;
        vector_normalize(state.dir);
    }

    //! Change the direction through an angle
    void change_direction(
            double       costheta,
            double       phi,
            Geo_State_t& state)
    {
        cartesian_vector_transform(costheta, phi, state.dir);
    }

    //! Reflect the direction at a reflecting surface.
    bool reflect(Geo_State_t& state)
    {
        NOT_IMPLEMENTED("reflect in mesh_geometry");
    }

    //! Return the outward normal at the location dictated by the state.
    Space_Vector normal(const Geo_State_t& state) const
    {
        NOT_IMPLEMENTED("normal in mesh_geometry");
    }

    // >>> PUBLIC INTERFACE

    // Get the volumes.
    SP_Vec_Dbl get_cell_volumes();

    // If the particle is outside the geometry, find distance
    double distance_to_interior(Geo_State_t &state);

    //! Access the underlying mesh directly
    const Cartesian_Mesh& mesh() const { return d_mesh; }

  private:
    // >>> IMPLEMENTATION

    // Update state tracking information
    void update_state(Geo_State_t &state) const;

    //! Move a particle a distance \e d in the current direction.
    void move(double dist, Geo_State_t &state) const
    {
        using profugus::soft_equiv;
        using def::X; using def::Y; using def::Z;

        REQUIRE(dist >= 0.0);
        REQUIRE(soft_equiv(vector_magnitude(state.dir), 1.0, 1.0e-6));

        // advance the particle (unrolled loop)
        state.pos[X] += dist * state.dir[X];
        state.pos[Y] += dist * state.dir[Y];
        state.pos[Z] += dist * state.dir[Z];
    }
};

} // end namespace profugus

#endif // geometry_Mesh_Geometry_hh

//---------------------------------------------------------------------------//
//                 end of Mesh_Geometry.hh
//---------------------------------------------------------------------------//
