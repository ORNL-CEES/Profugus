//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/RTK_Geometry.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 25 10:02:33 2011
 * \brief  RTK_Geometry class definition.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_RTK_Geometry_hh
#define geometry_RTK_Geometry_hh

#include <cmath>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "utils/Constants.hh"
#include "utils/SP.hh"
#include "utils/Definitions.hh"
#include "utils/Vector_Functions.hh"
#include "geometry/Definitions.hh"
#include "RTK_State.hh"
#include "RTK_Cell.hh"
#include "RTK_Array.hh"
#include "geometry/MOC_Geometry.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class RTK_Geometry
 * \brief Defines geometry implementations for the RTK MC geometry.
 *
 * The profugus::RTK_Array and profugus::RTK_Cell classes allow clients to build
 * hierarchical arrays of pin-cells, arrays of arrays of pin-cells, etc.  In
 * this way, simple LWR reactor geometries can be constructed quickly.
 * However, in order to use these classes in the profugus MC framework, they
 * must be defined as a profugus::Geometry.  The RTK_Geometry class provides the
 * correct publicly derived interface to profugus::Geometry so that the profugus
 * MC classes can use this as a geometry implementation.  This class is
 * parameterized on Array so that different types of geometries can be
 * constructed from the same code, ie.
 * \code
 typedef RTK_Geometry< RTK_Array<RTK_Cell> >              RTK_Lattice;
 typedef RTK_Geometry< RTK_Array< RTK_Array<RTK_Cell> > > RTK_Core;

 // make an RTK_Core core reactor geometry
 typedef profugus::RTK_Core        Core_Geometry;
 typedef Core_Geometry::Base     Geometry;
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
 * \sa profugus::RTK_Array, profugus::RTK_Cell, profugus::Geometry
 */
/*!
 * \example geometry/rtk/test/tstRTK_Lattice.cc
 * \example geometry/rtk/test/tstRTK_Core.cc
 *
 * Test of RTK_Geometry implementations.
 */
//===========================================================================//

template<class Array>
class RTK_Geometry : public MOC_Geometry<RTK_State>
{
    typedef MOC_Geometry<RTK_State>     Base;
  public:
    //@{
    //! Typedefs.
    typedef Array                       Array_t;
    typedef nemesis::SP<Array_t>        SP_Array;
    typedef RTK_State                   Geo_State_t;
    typedef typename Base::Space_Vector Space_Vector;
    typedef typename Base::RCF_Vec_Int  RCF_Vec_Int;
    typedef typename Base::RCF_Vec_Dbl  RCF_Vec_Dbl;
    typedef def::Vec_Dbl                Vec_Dbl;
    typedef def::Vec_Int                Vec_Int;
    //@}

  private:
    // >>> DATA

    // Underlying object.
    SP_Array d_array;

    // Volumes of each cell
    Vec_Dbl d_volumes;

    // Vector of symmetry cell ids
    Vec_Int d_mapped_cells;

    // Level of array.
    const int d_level;

    // Modular interface is determined from database
    bool d_modular;

  public:
    // Constructor.
    explicit RTK_Geometry(SP_Array array);

    // Constructor with symmetry
    RTK_Geometry(SP_Array array, bool sym);

    // >>> DERIVED INTERFACE from Geometry

    //! Initialize a track.
    void initialize(const Space_Vector &r, const Space_Vector &direction,
                    Geo_State_t &state);

    //! Get distance to next boundary.
    double distance_to_boundary(Geo_State_t &state)
    {
        Require (d_array);
        d_array->distance_to_boundary(state.d_r, state.d_dir, state);
        return state.dist_to_next_region;
    }

    //! Move to and cross a cell surface (do not reflect the particle, but
    //! indicate that the particle is on a reflecting surface).
    void move_to_surface(Geo_State_t &state)
    {
        Require (d_array);

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
        Require (d_array);

        // move the particle
        move(d, state);

        // update the array state to clear any surface tags
        d_array->update_state(state);
    }

    //! Return the current cell ID
    cell_type cell(const Geo_State_t &state) const
    {
        // Note: large arrays may overflow signed int
        Ensure(d_array->cellid(state) >= 0);
        return d_array->cellid(state);
    }

    //! Return the current material ID
    matid_type matid(const Geo_State_t &state) const
    {
        // we need a better mechanism later on....TME
        return d_array->matid(state);
    }

    //! Return the state with respect to outer geometry boundary
    Boundary_State boundary_state(const Geo_State_t &state) const
    {
        if (state.escaping_face != Geo_State_t::NONE)
        {
            // if the particle has escaped indicate that the particle
            // is outside the geometry
            Check (state.exiting_level[d_level]
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

    // Set volumes by global id
    void set_volume(int global_cell_id, double volume)
    {
        Check (global_cell_id < d_volumes.size() );
        d_volumes[global_cell_id] = volume;
    }

    // >>> DERIVED INTERFACE from MOC_Geometry

    //! Return Block/Plane id.  If non-modular will return 0.
    size_type block(const Geo_State_t &state) const
    {
        if (!d_modular)
            return 0;
        else
            return block_cmfd( state );
    }

    //! Return modular block/plane id
    size_type block_cmfd(const Geo_State_t &state) const
    {
        return d_array->index(
            state.level_coord[ d_array->level() ][def::X],
            state.level_coord[ d_array->level() ][def::Y], 0 );
    }

    //! Number of cells
    long_type num_cells() const
    {
        return d_array->num_cells();
    }

    //! Get user-supplied cell "labels" (unavailable, unassigned RCF)
    RCF_Vec_Int get_cell_labels() const
    {
        return RCF_Vec_Int();
    }

    //! Get cell volumes (unavailable, unassigned RCF)
    RCF_Vec_Dbl get_cell_volumes() const
    {
        return RCF_Vec_Dbl();
    }

    //! Number of cells in a block
    long_type num_cells_block(size_type block_id) const
    {
        if ( !d_modular )
            return d_array->num_cells();
        else
            return num_cells_block_cmfd( block_id );
    }

    //! Number of cells in a modular block
    long_type num_cells_block_cmfd(size_type block_id) const
    {
        int Nx = d_array->size( def::X);
        int i = block_id % Nx;
        return d_array->object( i, (block_id - i) / Nx, 0).num_cells();
    }

    //! Number of blocks
    long_type num_blocks() const
    {
        if ( !d_modular)
            return 1;
        else
            return num_blocks_cmfd();
    }

    //! Number of modular blocks
    long_type num_blocks_cmfd() const
    {
        return d_array->size();
    }

    //! Number of blocks in a given direction
    long_type num_blocks(size_type xyz) const
    {
        if ( !d_modular)
            return 1;
        else
            return num_blocks_cmfd( xyz );
    }

    //! Number of modular blocks in a given direction
    long_type num_blocks_cmfd(size_type xyz) const
    {
        return d_array->size( xyz );
    }

    //! Volume of global cell  // NEED TO UPDATE
    // THIS IS 2D VOLUME( AREA)
    double volume(size_type global_cell_id) const
    {
        Check (global_cell_id < d_volumes.size());
        return d_volumes[global_cell_id];
    }

    // Return all volumes
    const Vec_Dbl & volumes() const
    {
        return d_volumes;
    }

    //! Dimensions of block (x or y)
    double block_delta(size_type xyz, size_type block_id) const
    {
        if ( !d_modular)
        {
            Check( block_id == 0 );
            return d_array->pitch(xyz);
        }
        else
            return block_delta_cmfd( xyz, block_id );
    }

    //! Dimension of modular block (x or y)
    double block_delta_cmfd(size_type xyz, size_type block_id) const
    {
        int Nx = d_array->size(def::X);
        int i = block_id % Nx;
        return d_array->object( i, (block_id - i) / Nx, 0).pitch(xyz);
    }

    //! True if the particle/ray will change units
    bool change_planes(const Geo_State_t &state) const
    {
        // Plane boundaries begin at id 1000 (internal boundaries are < 1000)
        return state.exiting_face >= 1000;
    }

    //! Return bottom corner of a given block
    Space_Vector corner() const
    {
        return d_array->corner();
    }

    // Return the outward normal.
    Space_Vector normal(const Geo_State_t &state) const;

    //! Pickle the state.
    /// The RTK state is always persistent, so this is a null-op.
    void pickle(Geo_State_t &state) { /*...*/ }

    //! Restore the geometry from a persistent state.
    /// All of the geometric state is stored in the state, so this is a
    /// null-op for this geometry.
    void restore(Geo_State_t &pickled_state) { /*...*/ }

    //! Set modular interface
    void set_modular(bool mod) { d_modular = mod; }

    // Sets symmetric/mapped cells for eighth core symmetry
    void set_mapped_cells();

    // >>> ACCESSORS

    //! Return the underlying array representation of objects.
    const Array_t& array() const { Require (d_array); return *d_array; }

    //! Returns mapped cells
    int mapped_cell( int index )
    {
        Require( index < d_mapped_cells.size() );
        return d_mapped_cells[index];
    }

    //! Return mapped cell vector
    const Vec_Int& mapped_cells()
    {
        Require( !d_mapped_cells.empty() );
        return d_mapped_cells;
    }

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
        Require (d >= 0.0);
        Require (nemesis::soft_equiv(
                     vector_magnitude(state.d_dir), 1.0, 1.0e-6));

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
typedef RTK_Geometry< RTK_Array<RTK_Cell> >              RTK_Lattice;
typedef RTK_Geometry< RTK_Array< RTK_Array<RTK_Cell> > > RTK_Core;
//@}

} // end namespace profugus

#endif // geometry_RTK_Geometry_hh

//---------------------------------------------------------------------------//
//              end of geometry/RTK_Geometry.hh
//---------------------------------------------------------------------------//
