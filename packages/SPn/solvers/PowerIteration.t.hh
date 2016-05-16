//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/PowerIteration.t.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  PowerIteration template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_solvers_PowerIteration_t_hh
#define SPn_solvers_PowerIteration_t_hh

#include <vector>

#include "comm/P_Stream.hh"
#include "PowerIteration.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Build a PowerIteration solver.
 */
template <class T>
PowerIteration<T>::PowerIteration( RCP_ParameterList db )
    : EigenvalueSolver<T>(db)
{
    b_label = "Power Iteration";
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Solver a linear system using bicgstab.
 */
template <class T>
void PowerIteration<T>::solve( double           &lambda,
                                   Teuchos::RCP<MV>  x )
{
    REQUIRE( !b_A.is_null() );
    REQUIRE( !x.is_null() );

    // Allocate necessary vectors
    Teuchos::RCP<MV> Ax = MVT::Clone(*x,1);
    Teuchos::RCP<MV> r  = MVT::Clone(*x,1);

    std::vector<double> tmp_dot(1); // Temp storage for dot product.
    std::vector<double> tmp_nrm(1); // Temp storage for vector norm.
    double res_norm;

    // Normalize initial vector, if it is all zeros set it to a constant
    MVT::MvNorm(*x,tmp_nrm);
    if( tmp_nrm[0] == SCT::zero() )
    {
        MVT::MvInit(*x,1.0);
        MVT::MvNorm(*x,tmp_nrm);
    }
    MVT::MvScale(*x,1.0/tmp_nrm[0]);

    b_converged = false;
    b_num_iters = 0;

    while( true )
    {
        b_num_iters++;

        // Ax = A*x
        OPT::Apply(*b_A,*x,*Ax);

        // Compute new eigenvalue as Rayleigh quotient
        MVT::MvDot(*x,*Ax,tmp_dot);
        lambda = tmp_dot[0];

        // Compute residual
        MVT::MvAddMv(1.0,*Ax,-lambda,*x,*r);
        MVT::MvNorm(*r,tmp_nrm);
        res_norm = tmp_nrm[0];

        if( b_verbosity >= EigenvalueSolver<T>::LOW )
        {
            profugus::pout << " Power Iteration eigenvalue at iteration "
                           << b_num_iters << " is "
                           << profugus::fixed << profugus::setprecision(8)
                           << lambda << " with a residual norm of "
                           << profugus::scientific << profugus::setprecision(3)
                           << res_norm << profugus::endl;
        }

        // Set new eigenvector and normalize
        MVT::Assign(*Ax,*x);
        MVT::MvNorm(*x,tmp_nrm);
        MVT::MvScale(*x,1.0/tmp_nrm[0]);

        // Check for convergence
        if( res_norm < b_tolerance )
        {
            b_converged = true;
            break;
        }

        if( b_num_iters >= b_max_iters )
            break;
    }

    if( b_verbosity >= EigenvalueSolver<T>::MEDIUM )
    {
        profugus::pout << "+++ Power Iteration converged in "
                       << b_num_iters << " iterations." << profugus::endl;
    }
}

} // end namespace profugus

#endif // SPn_solvers_PowerIteration_t_hh

//---------------------------------------------------------------------------//
//                 end of PowerIteration.t.hh
//---------------------------------------------------------------------------//
