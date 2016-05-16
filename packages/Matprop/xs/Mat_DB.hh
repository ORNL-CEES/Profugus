//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Matprop/xs/Mat_DB.hh
 * \author Thomas M. Evans
 * \date   Mon Feb 10 19:50:43 2014
 * \brief  Mat_DB class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Matprop_xs_Mat_DB_hh
#define Matprop_xs_Mat_DB_hh

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "Teuchos_RCP.hpp"

#include "XS.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Mat_DB
 * \brief Material database.
 *
 * This class binds a cross section database to a list of cells.
 */
/*!
 * \example xs/test/tstMat_DB.cc
 *
 * Test of Mat_DB.
 */
//===========================================================================//

class Mat_DB
{
  public:
    //@{
    //! Typedefs.
    typedef XS                 XS_t;
    typedef Teuchos::RCP<XS_t> RCP_XS;
    typedef def::Vec_Int       Vec_Int;
    //@}

  private:
    // >>> DATA

    // Cells list of material ids.
    Vec_Int d_matids;

    // Cross sections.
    RCP_XS d_xs;

  public:
    // Constructor.
    Mat_DB();

    // >>> PUBLIC INTERFACE

    // Set cross section database and number of cells.
    void set(RCP_XS xs, int nc = 0);

    //! Mutable access to the material id.
    int& matid(int cell)
    {
        REQUIRE(cell < d_matids.size());
        return d_matids[cell];
    }

    //! Constant access to the material id.
    int matid(int cell) const
    {
        REQUIRE(cell < d_matids.size());
        return d_matids[cell];
    }

    // Assignment from a vector of matids.
    void assign(const Vec_Int &matids);

    // >>> ACCESSORS

    //! Get the cross section database.
    const XS_t& xs() const { REQUIRE(!d_xs.is_null()); return *d_xs; }

    //! Get the matids.
    const Vec_Int& matids() const { return d_matids; }

    //! Number of cells (zones) in the database.
    int num_cells() const { return d_matids.size(); }
};

} // end namespace profugus

#endif // Matprop_xs_Mat_DB_hh

//---------------------------------------------------------------------------//
//                 end of Mat_DB.hh
//---------------------------------------------------------------------------//
