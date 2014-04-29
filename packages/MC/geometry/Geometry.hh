//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   geometry/Geometry.hh
 * \author Thomas M. Evans
 * \date   Thu Dec 09 10:57:45 2010
 * \brief  Geometry class definition.
 * \note   Copyright (C) 2010 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef geometry_Geometry_hh
#define geometry_Geometry_hh

#include <vector>

#include "utils/Definitions.hh"
#include "utils/RCF.hh"
#include "Definitions.hh"

namespace denovo
{

//===========================================================================//
/*!
 * \class Geometry
 * \brief Defines interface for MC geometry packages.
 */
//===========================================================================//

template<class T>
class Geometry
{
  public:
    //@{
    //! Typedefs.
    typedef T                        Geo_State_t;
    typedef def::Space_Vector        Space_Vector;
    typedef geometry::Boundary_State Boundary_State;
    typedef geometry::matid_type     matid_type;
    typedef geometry::cell_type      cell_type;

    typedef nemesis::RCF<std::vector<int> >    RCF_Vec_Int; // cell labels
    typedef nemesis::RCF<std::vector<double> > RCF_Vec_Dbl;
    //@}

  public:
    // Virtual destructor.
    virtual ~Geometry() { /* * */ }

    // >>> PUBLIC INTERFACE

    //! Initialize track.
    virtual void initialize(const Space_Vector &r,
                            const Space_Vector &direction,
                            Geo_State_t        &state) = 0;

    //! Get distance to next boundary.
    virtual double distance_to_boundary(Geo_State_t &state) = 0;

    //! Move to and cross a surface in the current direction.
    virtual void move_to_surface(Geo_State_t &state) = 0;

    //! Move an arbitrary distance \e d to a point in the current direction.
    virtual void move_to_point(double d, Geo_State_t &state) = 0;

    //! Return the current material ID
    virtual matid_type matid(const Geo_State_t &state) const = 0;

    //! Return the current cell ID
    virtual cell_type cell(const Geo_State_t &state) const = 0;

    //! Return the state with respect to outer geometry boundary
    virtual Boundary_State boundary_state(const Geo_State_t &state) const = 0;

    //! Return the current position.
    virtual Space_Vector position(const Geo_State_t &state) const = 0;

    //! Return the current direction.
    virtual Space_Vector direction(const Geo_State_t &state) const = 0;

    //! Change the direction to \p new_direction.
    virtual void change_direction(const Space_Vector &new_direction,
                                  Geo_State_t &state) = 0;

    //! Change the direction through
    virtual void change_direction(double costheta, double phi,
                                  Geo_State_t &state) = 0;

    //! Reflect the direction at a reflecting surface.
    virtual bool reflect(Geo_State_t &state) = 0;

    //! Return the outward normal at the location dictated by the state.
    virtual Space_Vector normal(const Geo_State_t &state) const = 0;

    //! Pickle the state to make it persistent.
    /// Store the minimum amount of information in the state to reinitialize
    /// the geometry at a later time.
    virtual void pickle(Geo_State_t &state) = 0;

    //! Restore the geometry from a persistent ("pickled") state.
    virtual void restore(Geo_State_t &pickled_state) = 0;

    //! Get user-supplied cell "labels" (integers for now)
    virtual RCF_Vec_Int get_cell_labels() const = 0;

    //! Calculate cell volumes (potentially expensive)
    virtual RCF_Vec_Dbl get_cell_volumes() const = 0;
};

} // end namespace denovo

#endif // geometry_Geometry_hh

//---------------------------------------------------------------------------//
//              end of geometry/Geometry.hh
//---------------------------------------------------------------------------//
