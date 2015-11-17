//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/Geometry.hh
 * \author Thomas M. Evans
 * \date   Tuesday April 29 16:43:25 2014
 * \brief  Geometry class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_Geometry_hh
#define geometry_Geometry_hh

#include <cmath>
#include <memory>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "utils/Constants.hh"
#include "utils/Vector_Functions.hh"
#include "RTK_State.hh"
#include "RTK_Cell.hh"
#include "RTK_Array.hh"
#include "Tracking_Geometry.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Geometry
 * \brief Defines geometry implementations for the RTK MC geometry.
 *
 * The profugus::RTK_Array and profugus::RTK_Cell classes allow clients to
 * build hierarchical arrays of pin-cells, arrays of arrays of pin-cells, etc.
 * In this way, simple LWR reactor geometries can be constructed quickly.
 * However, in order to use these classes in the profugus MC framework, they
 * must be defined as a profugus::Geometry.  The Geometry class provides the
 * correct publicly derived interface to profugus::Geometry so that the
 * profugus MC classes can use this as a geometry implementation.  This class
 * is parameterized on Array so that different types of geometries can be
 * constructed from the same code, ie.
 * \code
 typedef Geometry< RTK_Array<RTK_Cell> >              Lattice;
 typedef Geometry< RTK_Array< RTK_Array<RTK_Cell> > > Core;

 // make an RTK_Core core reactor geometry
 typedef profugus::Core          Core_Geometry;
 typedef Core_Geometry::Array_t  Core_t;
 typedef Core_t::Object_t        Lattice_t;
 typedef Lattice_t::Object_t     Pin_Cell_t;
 typedef Core_Geometry::SP_Array SP_Core;
 typedef Core_t::SP_Object       SP_Lattice;
 typedef Lattice_t::SP_Object    SP_Pin_Cell;

 SP_Core core(new Core_t(...));
 // ...

 Core_Geometry rtk_core(core);
 * \endcode
 * See the tests for more examples.
 *
 * \sa profugus::RTK_Array, profugus::RTK_Cell
 */
/*!
 * \example geometry/rtk/test/tstLattice.cc
 * \example geometry/rtk/test/tstCore.cc
 *
 * Test of Geometry implementations.
 */
//===========================================================================//

template<class Array>
class Geometry : public Tracking_Geometry<RTK_State>
{
  public:
    //@{
    //! Typedefs.
    typedef Array                    Array_t;
    typedef std::shared_ptr<Array_t> SP_Array;
    typedef def::Vec_Dbl             Vec_Dbl;
    typedef def::Vec_Int             Vec_Int;
    //@}

  private:
    // >>> DATA

    // Underlying object.
    SP_Array d_array;

    // Volumes of each cell
    Vec_Dbl d_volumes;

    // Level of array.
    const int d_level;

  public:
    // Constructor.
    explicit Geometry(SP_Array array);

    // >>> DERIVED INTERFACE from Geometry_Base

    //! Initialize a track.
    void initialize(const Space_Vector &r, const Space_Vector &direction,
                    Geo_State_t &state) const;

    //! Get distance to next boundary.
    double distance_to_boundary(Geo_State_t &state)
    {
        REQUIRE(d_array);
        d_array->distance_to_boundary(state.d_r, state.d_dir, state);
        return state.dist_to_next_region;
    }

    //! Move to and cross a cell surface (do not reflect the particle, but
    //! indicate that the particle is on a reflecting surface).
    void move_to_surface(Geo_State_t &state)
    {
        REQUIRE(d_array);

        // move the particle
        move(state.dist_to_next_region, state);

        // process the particle through the surface
        d_array->cross_surface(state.d_r, state);
    }

    //! Move the particle to a point in the current direction.
    /// Clear any surface tags; the final point should \b not be a boundary
    /// surface.
    void move_to_point(double d, Geo_State_t &state)
    {
        REQUIRE(d_array);

        // move the particle
        move(d, state);

        // update the array state to clear any surface tags
        d_array->update_state(state);
    }

    //! Number of cells (excluding "outside" cell)
    geometry::cell_type num_cells() const { return d_array->num_cells(); }

    // Get the volume for a cell
    double cell_volume(int cellid) const
    {
        REQUIRE(cellid < num_cells());
        CHECK(cellid < d_volumes.size());
        double vol = d_volumes[cellid];
        ENSURE(vol >= 0.0);
        return vol;
    }

    //! Return the current cell ID
    geometry::cell_type cell(const Geo_State_t &state) const
    {
        // Note: large arrays may overflow signed int
        ENSURE(d_array->cellid(state) >= 0);
        return d_array->cellid(state);
    }

    //! Return the cell ID from the given location
    geometry::cell_type cell(const Space_Vector &r) const
    {
        return Tracking_Geometry<RTK_State>::cell(r);
    }

    //! Return the current material ID
    geometry::matid_type matid(const Geo_State_t &state) const
    {
        // we need a better mechanism later on....TME
        return d_array->matid(state);
    }

    //! Return the material ID from the given location
    geometry::matid_type matid(const Space_Vector &r) const
    {
        return Tracking_Geometry<RTK_State>::matid(r);
    }

    //! Return the state with respect to outer geometry boundary
    geometry::Boundary_State boundary_state(const Geo_State_t &state) const
    {
        if (state.escaping_face != Geo_State_t::NONE)
        {
            // if the particle has escaped indicate that the particle
            // is outside the geometry
            CHECK(state.exiting_level[d_level]
                    || state.next_face == Geo_State_t::NONE);
            return geometry::OUTSIDE;
        }
        else if (state.reflecting_face != Geo_State_t::NONE)
        {
            // test of reflection on the given face
            return geometry::REFLECT;
        }
        return geometry::INSIDE;
    }

    //! Return the boundary state for a given position
    geometry::Boundary_State boundary_state(const Space_Vector &r) const
    {
        return Tracking_Geometry<RTK_State>::boundary_state(r);
    }

    //! Return the current position.
    Space_Vector position(const Geo_State_t &state) const { return state.d_r; }

    //! Return the current direction.
    Space_Vector direction(const Geo_State_t &state) const {return state.d_dir;}

    //! Change the particle direction.
    void change_direction(const Space_Vector &new_direction, Geo_State_t &state)
    {
        // update the direction
        state.d_dir = new_direction;

        // normalize the direction
        vector_normalize(state.d_dir);
    }

    // Change the direction through angles \f$(\theta,\phi)\f$.
    void change_direction(double costheta, double phi, Geo_State_t &state)
    {
        cartesian_vector_transform(costheta, phi, state.d_dir);
    }

    // Reflect the direction at a reflecting surface.
    bool reflect(Geo_State_t &state);

    // Return the outward normal.
    Space_Vector normal(const Geo_State_t &state) const;

    // >>> ACCESSORS

    //! Return the underlying array representation of objects.
    const Array_t& array() const { REQUIRE(d_array); return *d_array; }

  private:
    // >>> IMPLEMENTATION

    /*! \brief Move a particle a distance \e d in the current direction.
     *
     * A particle is moved through a distance \f$d\f$ according to,
     *   \f[
     *   \left(\begin{array}{l}
     *   x\\y\\z\end{array}\right) =
     *   \left(\begin{array}{l}
     *   x_o\\y_o\\z_o\end{array}\right) + d
     *   \left(\begin{array}{l}
     *   \Omega_x\\ \Omega_y\\ \Omega_z\end{array}\right)\:.
     *   \f]
     *
     * The direction vector of the particle must be a unit-vector, ie:
     *  \f$|\Omega| = 1\f$.
     */
    void move(double d, Geo_State_t &state)
    {
        REQUIRE(d >= 0.0);
        REQUIRE(soft_equiv(vector_magnitude(state.d_dir), 1.0, 1.0e-6));

        // advance the particle (unrolled loop)
        state.d_r[def::X] += d * state.d_dir[def::X];
        state.d_r[def::Y] += d * state.d_dir[def::Y];
        state.d_r[def::Z] += d * state.d_dir[def::Z];
    }

  private:
    // Lower and upper extents of the enclosed array
    Space_Vector d_lower;
    Space_Vector d_upper;
};

//---------------------------------------------------------------------------//
// GEOMETRY TYPES
//---------------------------------------------------------------------------//

//@{
//! Single-level lattice/core geometries.
typedef Geometry< RTK_Array<RTK_Cell> >              Lattice;
typedef Geometry< RTK_Array< RTK_Array<RTK_Cell> > > Core;
//@}

} // end namespace profugus

#endif // geometry_Geometry_hh

//---------------------------------------------------------------------------//
//              end of geometry/Geometry.hh
//---------------------------------------------------------------------------//
