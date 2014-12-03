//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BelosSolver.hh
 * \author Steven Hamilton
 * \brief  BelosSolver class declaration.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_BelosSolver_hh
#define solvers_BelosSolver_hh

#include <SPn/config.h>
#include "LinearSolver.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "BelosSolverManager.hpp"

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \class BelosSolver
 * \brief Use Belos to solve linear system.
 *
 * This class provides a thin wrapper around a Belos linear solver
 * (i.e. a Belos::SolverManager).  The "Belos" sublist of the provided
 * ParameterList is provided directly to the Belos::SolverFactory, so any
 * Belos solvers and corresponding options are automatically available.
 * Consult Belos documentation for full list of solver options.
 */
//---------------------------------------------------------------------------//
template <class T>
class BelosSolver : public LinearSolver<T>
{
  public:

    typedef typename T::MV MV;
    typedef typename T::OP OP;
    typedef typename T::ST ST;
    typedef typename T::LO LO;
    typedef typename T::GO GO;

    explicit BelosSolver(Teuchos::RCP<Teuchos::ParameterList> pl);

    // Set operator for linear system
    void set_operator(Teuchos::RCP<OP> A)
    {
        REQUIRE( A != Teuchos::null );
        b_A = A;
    }

    // Set preconditioner for linear solver
    void set_preconditioner(Teuchos::RCP<OP> P)
    {
        REQUIRE( P != Teuchos::null );
        d_P = P;
    }

    void solve(Teuchos::RCP<MV>       x,
               Teuchos::RCP<const MV> b);

  private:

    using LinearSolver<T>::b_db;
    using LinearSolver<T>::b_A;
    using LinearSolver<T>::b_tolerance;
    using LinearSolver<T>::b_max_iters;
    using LinearSolver<T>::b_num_iters;
    using LinearSolver<T>::b_converged;
    using LinearSolver<T>::b_label;
    using LinearSolver<T>::b_verbosity;

    // Preconditioner
    Teuchos::RCP<OP> d_P;

    // Belos solver
    Teuchos::RCP<Belos::SolverManager<ST,MV,OP> > d_belos_solver;
};

} // namespace profugus

#endif // MCREX_BelosSolver_hh

