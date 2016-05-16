//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/Dimensions.hh
 * \author Thomas M. Evans
 * \date   Tue Oct 23 21:17:05 2012
 * \brief  Dimensions class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_spn_Dimensions_hh
#define SPn_spn_Dimensions_hh

namespace profugus
{

//===========================================================================//
/*!
 * \class Dimensions
 * \brief Stores \f$SP_N\f$ solver dimensions.
 */
/*!
 * \example spn/test/tstDimensions.cc
 *
 * Test of Dimensions.
 */
//===========================================================================//

class Dimensions
{
  private:
    // >>> DATA

    // SPN order (1, 3, 5, 7).
    int d_spn;

    // Number of equations (1, 2, 3, 4).
    int d_num_eqs;

  public:
    // Constructor.
    Dimensions(int N);

    //! Number of SPN equations.
    int num_equations() const { return d_num_eqs; }

    //! Number of moments.
    int num_moments() const { return d_spn + 1; }

    // Maximum number of equations
    static int max_num_equations() { return 4; }
};

} // end namespace profugus

#endif // SPn_spn_Dimensions_hh

//---------------------------------------------------------------------------//
//              end of spn/Dimensions.hh
//---------------------------------------------------------------------------//
