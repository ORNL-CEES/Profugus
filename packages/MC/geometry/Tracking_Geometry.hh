//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/Tracking_Geometry.hh
 * \author Thomas M. Evans
 * \date   Tue Jul 22 13:03:26 2014
 * \brief  Tracking_Geometry class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_Tracking_Geometry_hh
#define geometry_Tracking_Geometry_hh

#include "utils/Definitions.hh"
#include "Definitions.hh"
#include "Bounding_Box.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Tracking_Geometry
 * \brief Base class for geometries that enable particle tracking.
 */
//===========================================================================//

template<class Geo_State>
class Tracking_Geometry
{
  public:
    //@{
    //! Class Typedefs.
    typedef Geo_State                       Geo_State_t;
    typedef def::Space_Vector               Space_Vector;
    typedef def::Vec_Dbl                    Vec_Dbl;
    //@}

  public:
    // Virtual destructor.
    virtual ~Tracking_Geometry() = 0;

    // >> INTERFACE

    //! Initialize track.
    virtual void initialize(const Space_Vector& r,
                            const Space_Vector& direction,
                            Geo_State_t       & state) const = 0;

    //! Get distance to next boundary.
    virtual double distance_to_boundary(Geo_State_t& state) = 0;

    //! Move to and cross a surface in the current direction.
    virtual void move_to_surface(Geo_State_t& state) = 0;

    //! Move a distance \e d to a point in the current direction.
    virtual void move_to_point(double d, Geo_State_t& state) = 0;

    //! Number of cells (excluding "outside" cell)
    virtual geometry::cell_type num_cells() const = 0;

    //! Return the current cell ID, valid only when inside the mesh
    virtual geometry::cell_type cell(const Geo_State_t &state) const = 0;

    //! Return the cell ID from a given position inside the mesh
    inline geometry::cell_type cell(const Space_Vector &r) const;

    //! Return the current material ID
    virtual geometry::matid_type matid(const Geo_State_t &state) const = 0;

    //! Return the mateiral ID from a given position inside the mesh
    inline geometry::matid_type matid(const Space_Vector &r) const;

    //! Return the state with respect to outer geometry boundary
    virtual geometry::Boundary_State boundary_state(const Geo_State_t &state)
        const = 0;

    //! Return the state with respect to outer geometry boundary for a point
    inline geometry::Boundary_State boundary_state(const Space_Vector &r) const;

    //! Return the current position.
    virtual Space_Vector position(const Geo_State_t& state) const = 0;


    //! Return the current direction.
    virtual Space_Vector direction(const Geo_State_t& state) const = 0;

    //! Change the direction to \p new_direction.
    virtual void change_direction(const Space_Vector& new_direction,
                                  Geo_State_t& state) = 0;

    //! Change the direction through an angle
    virtual void change_direction(double costheta, double phi,
                                  Geo_State_t& state) = 0;

    //! Reflect the direction at a reflecting surface.
    virtual bool reflect(Geo_State_t& state) = 0;

    //! Return the outward normal at the location dictated by the state.
    virtual Space_Vector normal(const Geo_State_t& state) const = 0;

    //! Get all cell volumes
    virtual const Vec_Dbl &cell_volumes() const = 0;

    //! Get bounding box
    virtual Bounding_Box get_extents() const = 0;

    //! Get bounding box for a cell
    virtual Bounding_Box get_cell_extents(geometry::cell_type cell) const = 0;
};

//---------------------------------------------------------------------------//
// TEMPLATE MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Geo_State>
Tracking_Geometry<Geo_State>::~Tracking_Geometry()
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the cell that contains a given position.
 */
template<class Geo_State>
geometry::cell_type
Tracking_Geometry<Geo_State>::cell(const Space_Vector &r) const
{
    // Create a dummy angle
    Space_Vector dummy_dir(1.0, 0.0, 0.0);

    // Create a temporary Geo state
    Geo_State_t geo_state;

    // Initialize the Geo state with the given position
    this->initialize(r, dummy_dir, geo_state);

    // Return the cell from the geometry state
    return this->cell(geo_state);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the matid for the given position.
 */
template<class Geo_State>
geometry::cell_type
Tracking_Geometry<Geo_State>::matid(const Space_Vector &r) const
{
    // Create a dummy angle
    Space_Vector dummy_dir(1.0, 0.0, 0.0);

    // Create a temporary Geo state
    Geo_State_t geo_state;

    // Initialize the Geo state with the given position
    this->initialize(r, dummy_dir, geo_state);

    // Return the cell from the geometry state
    return this->matid(geo_state);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the boundary state for a given position.
 */
template<class Geo_State>
geometry::Boundary_State
Tracking_Geometry<Geo_State>::boundary_state(const Space_Vector &r) const
{
    // Create a dummy angle
    Space_Vector dummy_dir(1.0, 0.0, 0.0);

    // Create a temporary Geo state
    Geo_State_t geo_state;

    // Initialize the Geo state with the given position
    this->initialize(r, dummy_dir, geo_state);

    // Return the cell from the geometry state
    return this->boundary_state(geo_state);
}

} // end namespace profugus

#endif // geometry_Tracking_Geometry_hh

//---------------------------------------------------------------------------//
//                 end of Tracking_Geometry.hh
//---------------------------------------------------------------------------//
