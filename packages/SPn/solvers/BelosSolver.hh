//----------------------------------*-C++-*----------------------------------//
/*!\file   BelosSolver.hh
 * \author Steven Hamilton
 * \brief  BelosSolver class declaration.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_solvers_BelosSolver_hh
#define SPn_solvers_BelosSolver_hh

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

    void set_tolerance(double tol)
    {
        REQUIRE( tol > 0.0 );
        b_tolerance = tol;
        Teuchos::sublist(b_db,"Belos")->set("Convergence Tolerance",tol);
    }

    void set_max_iters(int max_iters)
    {
        REQUIRE( max_iters > 0 );
        b_max_iters = max_iters;
        Teuchos::sublist(b_db,"Belos")->set("Maximum Iterations",max_iters);
    }

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

    // Type of solver (e.g. "Block GMRES")
    std::string d_belos_type;

    // List of Belos-specific entries
    Teuchos::RCP<Teuchos::ParameterList> d_belos_pl;

    // Belos solver
    Teuchos::RCP<Belos::SolverManager<ST,MV,OP> > d_belos_solver;
};

} // namespace profugus

#endif // SPn_solvers_BelosSolver_hh

