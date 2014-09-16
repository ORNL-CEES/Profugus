//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/RayleighQuotient.t.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Mon Feb 24 13:29:04 2014
 * \brief  RayleighQuotient member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_RayleighQuotient_t_hh
#define solvers_RayleighQuotient_t_hh

#include "harness/Warnings.hh"
#include "comm/P_Stream.hh"
#include "RayleighQuotient.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Build a RayleighQuotient solver.
 */
template <class T>
RayleighQuotient<T>::RayleighQuotient( RCP_ParameterList db )
    : Base(db)
{
    b_label           = "Rayleigh Quotient";
    d_use_fixed_shift = db->get("use_fixed_shift", false);

    if( d_use_fixed_shift )
    {
        INSIST( db->isParameter("eig_shift"),
                "Must specify 'eig_shift' if fixed shift is requested." );
        d_fixed_shift = 1.0/db->get<double>("eig_shift");
        VALIDATE( d_fixed_shift > 0.0,
                  "Value of 'eig_shift' must be positive." );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve an eigenvalue problem using Rayleigh quotient iteration.
 * Note that this solver is intended for solving the k-eigenvalue problem
 */
template <class T>
void RayleighQuotient<T>::solve( double           &keff,
                                     Teuchos::RCP<MV>  x )
{
    REQUIRE( !b_A.is_null() );
    REQUIRE( !d_Op.is_null() );
    REQUIRE( !x.is_null() );

    // Allocate necessary vectors
    Teuchos::RCP<MV> Ax = MVT::Clone(*x,1);
    Teuchos::RCP<MV> Bx = MVT::Clone(*x,1);
    Teuchos::RCP<MV> r  = MVT::Clone(*x,1);

    std::vector<double> tmp_dot(1); // Temp storage for dot product.
    std::vector<double> tmp_nrm(1); // Temp storage for vector norm.
    double res_norm, num, denom;

    bool have_B = !d_B.is_null();

    // Normalize initial vector, if it is all zeros set it to a constant
    MVT::MvNorm(*x,tmp_nrm);
    if( tmp_nrm[0] == SCT::zero() )
    {
        MVT::MvInit(*x,1.0);
        MVT::MvNorm(*x,tmp_nrm);
    }
    MVT::MvScale(*x,1.0/tmp_nrm[0]);

    // Initialize Bx
    if( have_B )
    {
        OPT::Apply(*d_B,*x,*Bx);
    }
    else
    {
        std::vector<int> ind(1);
        ind[0] = 0;
        MVT::SetBlock(*x,ind,*Bx);
    }

    double tol;
    double lambda = 1.0 / keff;
    if( d_use_fixed_shift )
        d_Op->set_shift(d_fixed_shift);

    b_converged = false;
    b_num_iters = 0;

    while( true )
    {
        b_num_iters++;

        // Set shift in inverse operator
        if( !d_use_fixed_shift )
            d_Op->set_shift(lambda);

        // x = (A-lambda*B)^{-1}*Bx
        OPT::Apply(*d_Op,*Bx,*x);

        // Normalize eigenvector
        MVT::MvNorm(*x,tmp_nrm);
        MVT::MvScale(*x,1.0/tmp_nrm[0]);

        // Apply operators
        OPT::Apply(*b_A,*x,*Ax);
        if( have_B )
        {
            OPT::Apply(*d_B,*x,*Bx);
        }
        else
        {
            std::vector<int> ind(1);
            ind[0] = 0;
            MVT::SetBlock(*x,ind,*Bx);
        }

        // Compute dot products
        MVT::MvDot(*x,*Bx,tmp_dot);
        num = tmp_dot[0];
        MVT::MvDot(*x,*Ax,tmp_dot);
        denom = tmp_dot[0];

        // New eigenvalue: <x,Bx>/<x,Ax>
        keff= num / denom;
        lambda = denom / num;

        // Compute residual
        MVT::MvAddMv(1.0,*Ax,-1.0/keff,*Bx,*r);
        MVT::MvNorm(*r,tmp_nrm);
        res_norm = tmp_nrm[0];

        if( b_verbosity >= EigenvalueSolver<T>::LOW )
        {
            profugus::pout << " Rayleigh Quotient eigenvalue at iteration "
                           << profugus::fixed << profugus::setprecision(8)
                           << b_num_iters << " is " << keff
                           << " with a residual norm of "
                           << profugus::scientific << profugus::setprecision(4)
                           << res_norm << profugus::endl;
        }

        // Check for convergence
        if( res_norm < b_tolerance )
        {
            b_converged = true;
            break;
        }

        if( b_num_iters >= b_max_iters )
            break;

    }

    if( b_verbosity >= EigenvalueSolver<T>::LOW )
    {
        profugus::pout << "+++ Rayleigh Quotient Iteration converged in "
                       << b_num_iters << " iterations." << profugus::endl;
    }
}

} // end namespace profugus

#endif // solvers_RayleighQuotient_t_hh

//---------------------------------------------------------------------------//
//                 end of RayleighQuotient.cc
//---------------------------------------------------------------------------//
