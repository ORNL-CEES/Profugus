//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/PowerIteration.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  PowerIteration class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_PowerIteration_hh
#define solvers_PowerIteration_hh

#include "EigenvalueSolver.hh"

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"

namespace profugus
{

//===========================================================================//
/*!
 * \class PowerIteration
 * \brief PowerIteration solver.
 *
 * The constructor takes in an SP to a Std_DB.  The following entries are
 * significant:
 *  - ``tolerance'' The tolerance for the relative residual.
 *  - ``max_itr''   The maximum number of iterations to be performed.
 *
 * \sa PowerIteration.t.hh for detailed descriptions.
 */
/*!
 * \example solvers/test/tstPowerIteration.cc
 *
 * Test of PowerIteration.
 */
//===========================================================================//

template <class MV, class OP>
class PowerIteration : public EigenvalueSolver<MV,OP>
{
  public:
    //@{
    //! Typedefs.
    typedef EigenvalueSolver<MV,OP>               Base;
    typedef typename Base::ParameterList          ParameterList;
    typedef typename Base::RCP_ParameterList      RCP_ParameterList;
    typedef Teuchos::ScalarTraits<double>         SCT;
    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    typedef Anasazi::OperatorTraits<double,MV,OP> OPT;
    //@}

    // Constructor
    PowerIteration( RCP_ParameterList db );

    // Solve
    void solve( double           &lambda,
                Teuchos::RCP<MV>  x );

  private:

    using EigenvalueSolver<MV,OP>::b_db;
    using EigenvalueSolver<MV,OP>::b_A;
    using EigenvalueSolver<MV,OP>::b_tolerance;
    using EigenvalueSolver<MV,OP>::b_num_iters;
    using EigenvalueSolver<MV,OP>::b_max_iters;
    using EigenvalueSolver<MV,OP>::b_converged;
    using EigenvalueSolver<MV,OP>::b_label;
    using EigenvalueSolver<MV,OP>::b_verbosity;
};

} // end namespace profugus

#endif // solvers_PowerIteration_hh

//---------------------------------------------------------------------------//
//                 end of PowerIteration.hh
//---------------------------------------------------------------------------//
