//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AleaSolver.cc
 * \author Steven Hamilton
 * \brief  Perform Adjoint Monte Carlo on linear system
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "AleaSolver.hh"
#include "comm/global.hh"
#include "harness/DBC.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param A  Problem matrix
 * \param pl ParameterList
 *
 * Behavior is controlled by following PL entries:
 *  - max_iterations(int) : >0 (1000)
 *  - tolerance(double)   : >0.0 (1.0e-6)
 *  - verbosity(string)   : "none", ("low"), "medium", "high", "debug"
 */
//---------------------------------------------------------------------------//
AleaSolver::AleaSolver(Teuchos::RCP<const MATRIX> A,
                           Teuchos::RCP<Teuchos::ParameterList> pl )
  : b_A(A)
  , b_pl(pl)
{
    // Set default parameters on pl
    if( !pl->isType<int>("max_iterations") )
        pl->set<int>("max_iterations",1000);

    if( !pl->isType<MAGNITUDE>("tolerance") )
        pl->set<MAGNITUDE>("tolerance",1.0e-6);

    // Set the verbosity based on value in this PL
    if( !pl->isType<std::string>("verbosity") )
        pl->set<std::string>("verbosity","low");

    // Set class members
    setParameters(pl);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply Operator
 *
 *  This function partially implements the Tpetra Operator interface.
 *  In particular, it correctly deals with the scalar parameters alpha
 *  and beta, then passes the ``raw'' arguments x and y to applyImpl in the
 *  derived class.  This is done to avoid every derived class having the
 *  same set of checks for early exit and combining intermediate vectors at
 *  the end.  If a derived class can truly implement y = alpha*A*x + beta*y
 *  more efficiently than an operator apply followed by an axpy, it is free
 *  to overload this function directly.  Otherwise, derived classes need only
 *  implement the applyImpl function.
 */
//---------------------------------------------------------------------------//
void AleaSolver::apply(const MV &x, MV &y, Teuchos::ETransp mode,
                         SCALAR alpha, SCALAR beta) const
{
    REQUIRE( mode == Teuchos::NO_TRANS || this->hasTransposeApply() );
    REQUIRE( x.getNumVectors()  == y.getNumVectors() );
    REQUIRE( x.getLocalLength() == y.getLocalLength() );

    if( beta == SCALAR_TRAITS::zero() )
    {
        y.putScalar(0.0);
    }

    // Early exit if alpha is 0
    if( alpha == SCALAR_TRAITS::zero() )
    {
        y.scale(beta);
        return;
    }

    // Temporary vector to hold raw apply
    MV tmp_y(y.getMap(),y.getNumVectors());

    // "Raw" apply
    this->applyImpl(x,tmp_y);

    // y = alpha*tmp_y + beta*y
    y.update(alpha,tmp_y,beta);
}

//---------------------------------------------------------------------------//
// PROTECTED MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Update parameters from pl
//---------------------------------------------------------------------------//
void AleaSolver::setParameters( Teuchos::RCP<Teuchos::ParameterList> pl )
{
    if( pl->isType<int>("max_iterations") )
    {
        b_max_iterations = pl->get<int>("max_iterations");
    }

    if( pl->isType<MAGNITUDE>("tolerance") )
    {
        b_tolerance = pl->get<MAGNITUDE>("tolerance");
    }

    if( pl->isType<std::string>("verbosity") )
    {
        std::string verbosity = pl->get<std::string>("verbosity");
        VALIDATE(verbosity=="none"   || verbosity=="low"  ||
                 verbosity=="medium" || verbosity=="high" ||
                 verbosity=="debug",
                 "Invalid verbosity specified.");
        if( verbosity == "none")
        {
            b_verbosity = NONE;
        }
        else if( verbosity == "low" )
        {
            b_verbosity = LOW;
        }
        else if( verbosity == "medium" )
        {
            b_verbosity = MEDIUM;
        }
        else if( verbosity == "high" )
        {
            b_verbosity = HIGH;
        }
        else if( verbosity == "debug" )
        {
            b_verbosity = DEBUG;
        }
        //
        // Silence processors other than rank 0
        if( profugus::node() != 0 )
            b_verbosity = NONE;
    }
}

} // namespace alea

