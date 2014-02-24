//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/ShiftedInverseOperator.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 13:38:42 2014
 * \brief  ShiftedInverseOperator class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_ShiftedInverseOperator_hh
#define solvers_ShiftedInverseOperator_hh

#include <Teuchos_RCP.hpp>

#include "InverseOperator.hh"
#include "ShiftedOperator.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class ShiftedInverseOperator
 * \brief Perform shift-and-invert operation.
 *
 * This class performs the operation \f$ y = (A - \lambda I)^{-1}x \f$,
 * or \f$ y = (A - \lambda B)^{-1} Bx \f$.
 * \brief
 */
/*!
 * \example solvers/test/tstShiftedInverseOperator.cc
 *
 * Test of ShiftedInverseOperator.
 */
//===========================================================================//

class ShiftedInverseOperator : public InverseOperator
{
  private:
    typedef InverseOperator Base;

    // >>> DATA
    Teuchos::RCP<ShiftedOperator> d_operator;
    double d_shift;

  public:

    // Constructor.
    // Read Denovo database entries for solver parameters.
    explicit ShiftedInverseOperator(RCP_ParameterList db);

    void set_operator( RCP_Operator A);
    void set_rhs_operator( RCP_Operator A);

    // Set shift
    void set_shift( double shift )
    {
        d_operator->set_shift(shift);
    }

    // Apply (solve linear system)
    int Apply(const MV &x, MV &y) const;
};

} // end namespace profugus

#endif // solvers_ShiftedInverseOperator_hh

//---------------------------------------------------------------------------//
//                 end of ShiftedInverseOperator.hh
//---------------------------------------------------------------------------//
