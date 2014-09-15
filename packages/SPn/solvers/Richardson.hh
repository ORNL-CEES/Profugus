//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Richardson.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 12:06:57 2014
 * \brief  Richardson class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_Richardson_hh
#define solvers_Richardson_hh

#include "harness/DBC.hh"
#include "LinearSolver.hh"

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"

namespace profugus
{

//===========================================================================//
/*!
 * \class Richardson
 * \brief Richardson solver.
 *
 * The constructor takes in parameterlist.  The following entries are
 * significant:
 *  - ``tolerance'' The tolerance for the relative residual.
 *  - ``max_itr''   The maximum number of iterations to be performed.
 *  - ``damping''   Damping factor applied to iterates.
 * The memory requirement for this solver is three vectors (including the
 * solution vector and rhs, which are allocated outside of this class).
 * One additional vector is required if a preconditioner is used.
 *
 * \sa Richardson.t.hh for detailed descriptions.
 */
/*!
 * \example solvers/test/tstRichardson.cc
 *
 * Test of Richardson.
 */
//===========================================================================//

template <LinAlgType T>
class Richardson : public LinearSolver<T>
{
  public:
    //@{
    //! Typedefs.
    typedef typename LinAlgTypedefs<T>::MV        MV;
    typedef typename LinAlgTypedefs<T>::OP        OP;
    typedef LinearSolver<T>                       Base;
    typedef typename Base::ParameterList          ParameterList;
    typedef typename Base::RCP_ParameterList      RCP_ParameterList;
    typedef Teuchos::RCP<MV>                      RCP_MV;
    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    typedef Anasazi::OperatorTraits<double,MV,OP> OPT;
    //@}

  public:
    // Constructor
    Richardson( RCP_ParameterList db );

    // Solve
    void solve( Teuchos::RCP<MV>       x,
                Teuchos::RCP<const MV> b );

    //! Set the preconditioner.
    void set_preconditioner( Teuchos::RCP<OP> P )
    {
        REQUIRE( P != Teuchos::null );
        d_P = P;
    }

  private:

    Teuchos::RCP<OP> d_P;

    using LinearSolver<T>::b_db;
    using LinearSolver<T>::b_A;
    using LinearSolver<T>::b_tolerance;
    using LinearSolver<T>::b_num_iters;
    using LinearSolver<T>::b_max_iters;
    using LinearSolver<T>::b_converged;
    using LinearSolver<T>::b_label;
    using LinearSolver<T>::b_verbosity;

    double d_damping;
};

} // end namespace profugus

#endif // solvers_Richardson_hh

//---------------------------------------------------------------------------//
//                 end of Richardson.hh
//---------------------------------------------------------------------------//
