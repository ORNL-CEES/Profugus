//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Dimensions.cc
 * \author Thomas M. Evans
 * \date   Tue Oct 23 21:17:05 2012
 * \brief  Dimensions member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "harness/DBC.hh"
#include "Dimensions.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * The number of equations is equal to \f$(N+1)/2\f$.
 *
 * \param N SPN order (1, 3, 5, 7)
 */
Dimensions::Dimensions(int N)
    : d_spn(N)
{
    Require (N == 1 || N == 3 || N == 5 || N == 7);

    // set the number of equations
    d_num_eqs = (N + 1) / 2;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Dimensions.cc
//---------------------------------------------------------------------------//
