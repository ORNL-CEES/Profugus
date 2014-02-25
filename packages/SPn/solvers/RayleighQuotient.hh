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

class RayleighQuotient :
        public EigenvalueSolver<Epetra_MultiVector,Epetra_Operator>
{
    typedef EigenvalueSolver<Epetra_MultiVector,Epetra_Operator> Base;

  public:
    //@{
    //! Typedefs.
    typedef Epetra_MultiVector                    MV;
    typedef Epetra_Operator                       OP;
    typedef Teuchos::ScalarTraits<double>         SCT;
    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    typedef Anasazi::OperatorTraits<double,MV,OP> OPT;
    //@}

  public:
    // Constructor
    RayleighQuotient( RCP_ParameterList db );

    // Set shifted inverse operator
    void set_rhs_operator( Teuchos::RCP<OP> B )
    {
        Require( !B.is_null() );
        d_B = B;
    }

    // Set shifted inverse operator
    void set_shifted_operator( Teuchos::RCP<ShiftedInverseOperator> Op )
    {
        Require( !Op.is_null() );
        d_Op = Op;
    }

    // Solve
    void solve( double                &lambda,
                Teuchos::RCP<MV>       x );

  private:

    // Shifted inverse operator
    Teuchos::RCP<ShiftedInverseOperator> d_Op;
    Teuchos::RCP<OP>                     d_B;

    // Values for fixed (Wielandt) shift
    bool   d_use_fixed_shift;
    double d_fixed_shift;
};

} // end namespace profugus

#endif // solvers_RayleighQuotient_hh

//---------------------------------------------------------------------------//
//                 end of RayleighQuotient.hh
//---------------------------------------------------------------------------//
