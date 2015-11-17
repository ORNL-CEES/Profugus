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
    typedef Geo_State         Geo_State_t;
    typedef def::Space_Vector Space_Vector;
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

    //! Return the current material ID
    virtual geometry::matid_type matid(const Geo_State_t &state) const = 0;

    //! Return the state with respect to outer geometry boundary
    virtual geometry::Boundary_State boundary_state(const Geo_State_t &state)
        const = 0;

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

} // end namespace profugus

#endif // geometry_Tracking_Geometry_hh

//---------------------------------------------------------------------------//
//                 end of Tracking_Geometry.hh
//---------------------------------------------------------------------------//
