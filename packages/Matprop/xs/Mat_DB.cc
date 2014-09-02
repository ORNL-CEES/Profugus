//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/Mat_DB.cc
 * \author Thomas M. Evans
 * \date   Mon Feb 10 19:50:43 2014
 * \brief  Mat_DB member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "Mat_DB.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Mat_DB::Mat_DB()
{
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Set the number of cells and cross-section database.
 */
void Mat_DB::set(RCP_XS xs,
                 int    nc)
{
    REQUIRE(!xs.is_null());

    // cross section database
    d_xs = xs;

    // clear the existing database
    d_matids.clear();
    d_matids.resize(nc);
    std::fill(d_matids.begin(), d_matids.end(), -1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment from a vector of matids.
 *
 * \warning This will resize the internal number of cells to the input vector
 * size
 *
 * \param matids vector of material ids; the internal number of cells will be
 * set to the size of this vector
 */
void Mat_DB::assign(const Vec_Int &matids)
{
    REQUIRE(!matids.empty());
    d_matids = matids;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Mat_DB.cc
//---------------------------------------------------------------------------//
