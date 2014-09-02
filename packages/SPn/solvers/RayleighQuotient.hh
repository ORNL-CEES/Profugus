//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/RayleighQuotient.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Mon Feb 24 13:29:04 2014
 * \brief  RayleighQuotient class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_RayleighQuotient_hh
#define solvers_RayleighQuotient_hh

#include "EigenvalueSolver.hh"
#include "ShiftedInverseOperator.hh"

#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"

namespace profugus
{

//===========================================================================//
/*!
 * \class RayleighQuotient
 * \brief RayleighQuotient iteration solver.
 *
 * The constructor takes in an SP to a Std_DB.  The following entries are
 * significant:
 *  - ``tolerance'' The tolerance for the relative residual.
 *  - ``max_itr''   The maximum number of iterations to be performed.
 *
 * \sa RayleighQuotient.t.hh for detailed descriptions.
 */
/*!
 * \example solvers/test/tstRayleighQuotient.cc
 *
 * Test of RayleighQuotient.
 */
//===========================================================================//

template <class MV, class OP>
class RayleighQuotient : public EigenvalueSolver<MV,OP>
{
    typedef EigenvalueSolver<MV,OP> Base;

  public:
    //@{
    //! Typedefs.
    typedef Teuchos::ScalarTraits<double>         SCT;
    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    typedef Anasazi::OperatorTraits<double,MV,OP> OPT;
    typedef Teuchos::RCP<Teuchos::ParameterList>  RCP_ParameterList;
    //@}

  public:
    // Constructor
    RayleighQuotient( RCP_ParameterList db );

    // Set shifted inverse operator
    void set_rhs_operator( Teuchos::RCP<OP> B )
    {
        REQUIRE( !B.is_null() );
        d_B = B;
    }

    // Set shifted inverse operator
    void set_shifted_operator( Teuchos::RCP<ShiftedInverseOperator<MV,OP> > Op )
    {
        REQUIRE( !Op.is_null() );
        d_Op = Op;
    }

    // Solve
    void solve( double                &lambda,
                Teuchos::RCP<MV>       x );

  private:

    // Shifted inverse operator
    Teuchos::RCP<ShiftedInverseOperator<MV,OP> > d_Op;
    Teuchos::RCP<OP>                     d_B;

    // Values for fixed (Wielandt) shift
    bool   d_use_fixed_shift;
    double d_fixed_shift;

    using Base::b_A;
    using Base::b_converged;
    using Base::b_num_iters;
    using Base::b_verbosity;
    using Base::b_max_iters;
    using Base::b_tolerance;
    using Base::b_label;
};

} // end namespace profugus

#endif // solvers_RayleighQuotient_hh

//---------------------------------------------------------------------------//
//                 end of RayleighQuotient.hh
//---------------------------------------------------------------------------//
