//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/Richardson.t.hh
 * \author Thomas M. Evans
 * \date   Fri Feb 21 12:07:58 2014
 * \brief  Richardson template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_Richardson_t_hh
#define solvers_Richardson_t_hh

#include <vector>

#include "comm/P_Stream.hh"
#include "Richardson.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Build a native Richardson solver.
 */
template <LinAlgType T>
Richardson<T>::Richardson( RCP_ParameterList db )
    : LinearSolver<T>(db)
{
    d_damping = b_db->get("Damping Factor", 1.0);
    b_label   = "Profugus Richardson";
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Solve a linear system using damped Richardson iteration.
 */
template <LinAlgType T>
void Richardson<T>::solve( Teuchos::RCP<MV>       x,
                               Teuchos::RCP<const MV> b )
{
    REQUIRE( b_A != Teuchos::null );
    REQUIRE( x   != Teuchos::null );
    REQUIRE( b   != Teuchos::null );

    // Allocate necessary vectors
    Teuchos::RCP<MV> r = MVT::Clone(*x,1);

    Teuchos::RCP<MV> tmp;
    if( d_P != Teuchos::null )
    {
        tmp = MVT::Clone(*x,1);
    }

    std::vector<double> tmp_nrm(1); // Temp storage for vector norm.

    MVT::MvNorm(*b,tmp_nrm);
    double b_norm = tmp_nrm[0];
    double res_norm;

    b_converged = false;
    b_num_iters = 0;

    // If the RHS is zero, set solution to zero and return
    if( b_norm == 0.0 )
    {
        MVT::MvInit(*x,0.0);
        b_converged = true;
        return;
    }

    while( true )
    {
        // Compute residual
        if( d_P != Teuchos::null )
        {
            OPT::Apply(*b_A,*x,*tmp);
            MVT::MvAddMv(1.0,*b,-1.0,*tmp,*tmp);
            OPT::Apply(*d_P,*tmp,*r);
        }
        else
        {
            OPT::Apply(*b_A,*x,*r);
            MVT::MvAddMv(1.0,*b,-1.0,*r,*r);
        }

        // Check for convergence
        MVT::MvNorm(*r,tmp_nrm);
        res_norm = tmp_nrm[0];
        if( res_norm/b_norm < b_tolerance )
        {
            b_converged = true;
            break;
        }

        // Print status if requested
        if( b_verbosity >= LinearSolver<T>::MEDIUM )
        {
            profugus::pout << b_label << " residual norm at iteration "
                           << b_num_iters << " is "
                           << profugus::scientific << profugus::setprecision(3)
                           << res_norm/b_norm << profugus::endl;
        }

        // x = x + omega*r
        MVT::MvAddMv(1.0,*x,d_damping,*r,*x);

        b_num_iters++;

        // Check for max iterations
        if( b_num_iters >= b_max_iters )
            break;
    }

    // Print final status
    if( b_verbosity >= LinearSolver<T>::LOW )
    {
        if( b_converged )
        {
            profugus::pout << b_label << " converged after " << b_num_iters
                           << " iterations." << profugus::endl;
        }
        else
        {
            profugus::pout << b_label << " terminated after " << b_num_iters
                           << " iterations."  << profugus::endl
                           << " Final residual norm is " << res_norm/b_norm
                           << "." << profugus::endl
                           << " Requested tolerance is " << b_tolerance << "."
                           << profugus::endl;
        }
    }
}

} // end namespace profugus

#endif // solvers_Richardson_t_hh

//---------------------------------------------------------------------------//
//                 end of Richardson.t.hh
//---------------------------------------------------------------------------//
