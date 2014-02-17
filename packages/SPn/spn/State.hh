//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/State.hh
 * \author Thomas M. Evans
 * \date   Mon Feb 17 13:10:27 2014
 * \brief  State class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_State_hh
#define spn_State_hh

#include "Teuchos_ArrayView.hpp"
#include "Teuchos_RCP.hpp"

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "mesh/Mesh.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class State
 * \brief Holds unknowns (scalar flux) for SPn system/solver.
 */
/*!
 * \example spn/test/tstState.cc
 *
 * Test of State.
 */
//===========================================================================//

class State
{
  public:
    //@{
    //! Typedefs.
    typedef def::Vec_Dbl                     Field;
    typedef Teuchos::ArrayView<double>       View_Field;
    typedef Teuchos::ArrayView<const double> const_View_Field;
    typedef Teuchos::RCP<Mesh>               RCP_Mesh;
    //@}

  private:
    // >>> DATA

    // Mesh.
    RCP_Mesh d_mesh;

    // Number of groups.
    int d_num_groups;

    // Scalar fluxes -> 1-D array ordered index = n + Nc * g, where n is the
    // cell index, Nc is the number of cells, and g is the group index.
    Field d_flux;

  public:
    // Constructor.
    State(RCP_Mesh mesh, int num_groups = 1)
        : d_mesh(mesh)
        , d_num_groups(num_groups)
        , d_flux(mesh->num_cells() * d_num_groups, 0.0)
    {
        Require (!mesh.is_null());
        Require (d_num_groups > 0);
    }

    // >>> ACCESSORS

    //! Size of state.
    int size() const { return d_flux.size(); }

    //! Number of groups.
    int num_groups() const { return d_num_groups; }

    //! Number of (local) cells on the current domain.
    int num_cells() const { return d_mesh->num_cells(); }

    //! Get const-view to all fluxes.
    const_View_Field flux() const
    {
        return Teuchos::arrayViewFromVector(d_flux);
    }

    //! Get mutable-view to all fluxes.
    View_Field flux()
    {
        return Teuchos::arrayViewFromVector(d_flux);
    }

    //! Get a const-view to fluxes over a range of groups.
    const_View_Field flux(int first, int last) const
    {
        Require (1 + last - first > 0 && last - first < d_num_groups);
        return const_View_Field(&d_flux[0] + first * d_mesh->num_cells(),
                                d_mesh->num_cells() * (1 + last - first));
    }

    //! Get a mutable-view to fluxes over a range of groups.
    View_Field flux(int first, int last)
    {
        Require (1 + last - first > 0 && last - first < d_num_groups);
        return View_Field(&d_flux[0] + first * d_mesh->num_cells(),
                          d_mesh->num_cells() * (1 + last - first));
    }

    //! Return the index into the field.
    int index(int c, int g) const { return c + d_mesh->num_cells() * g; }
};

} // end namespace profugus

#endif // spn_State_hh

//---------------------------------------------------------------------------//
//                 end of State.hh
//---------------------------------------------------------------------------//
