//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/AndersonSolver.hh
 * \author Steven Hamilton
 * \date   Wed Apr 01 11:01:28 2015
 * \brief  AndersonSolver class definition.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_AndersonSolver_hh
#define solvers_AndersonSolver_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "NOX_Solver_Generic.H"

#include "harness/DBC.hh"
#include "ModelEvaluatorWrapper.hh"
#include "LinAlgTypedefs.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class AndersonSolver
 * \brief Solve nonlinear function using NOX Anderson Acceleration
 *
 *
 * \sa AndersonSolver.t.hh for detailed descriptions.
 */
/*!
 * \example solvers/test/tstAndersonSolver.cc
 *
 * Test of AndersonSolver.
 */
//===========================================================================//

template <class T>
class AndersonSolver
{
  public:

    typedef typename T::ST ST;
    typedef typename T::MV MV;
    typedef typename T::OP OP;

    AndersonSolver( Teuchos::RCP<const OP>               op,
                    Teuchos::RCP<Teuchos::ParameterList> pl );

    // Solve nonlinear problem
    void solve(Teuchos::RCP<MV> x);

  private:

    // Original user operator
    Teuchos::RCP<const OP> d_op;

    // Operator wrapped into Thyra ModelEvaluator interface
    Teuchos::RCP<ModelEvaluatorWrapper<T> > d_model_eval;

    // Solver parameters
    Teuchos::RCP<Teuchos::ParameterList> d_pl;

    // NOX Solver
    Teuchos::RCP<NOX::Solver::Generic> d_solver;
};

} // end namespace profugus

#endif // solvers_AndersonSolver_hh

//---------------------------------------------------------------------------//
//                 end of AndersonSolver.hh
//---------------------------------------------------------------------------//
