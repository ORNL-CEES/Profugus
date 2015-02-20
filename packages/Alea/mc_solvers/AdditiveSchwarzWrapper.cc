//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdditiveSchwarzWrapper.cc
 * \author Steven Hamilton
 * \brief  Wrap local solver into an AdditiveSchwarz solver.
 */
//---------------------------------------------------------------------------//

#include <iterator>
#include <string>

#include "AdditiveSchwarzWrapper.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param A  Problem matrix
 * \param pl ParameterList
 *
 * Behavior is controlled by the nested "AdditiveSchwarz" sublist.
 * Consult Ifpack2::AdditiveSchwarz documentation for options.
 */
//---------------------------------------------------------------------------//
AdditiveSchwarzWrapper::AdditiveSchwarzWrapper( Teuchos::RCP<const MATRIX> A,
        Teuchos::RCP<Ifpack2::Preconditioner<SCALAR,LO,GO,NODE> > prec,
        Teuchos::RCP<Teuchos::ParameterList> pl )
  : AleaSolver(A,pl)
{
    // Get belos pl
    Teuchos::RCP<Teuchos::ParameterList> as_pl =
        Teuchos::sublist(pl,"AdditiveSchwarz");

    // Override default parameters if present on sublist
    this->setParameters(as_pl);

    // Build AdditiveSchwarz
    d_schwarz = Teuchos::rcp( new Ifpack2::AdditiveSchwarz<CRS_MATRIX>(A) );
    d_schwarz->setParameterList(as_pl);
    d_schwarz->setInnerPreconditioner(prec);
    d_schwarz->compute();

    b_label = "AdditiveSchwarzWrapper";
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Solve linear system using Belos.
 */
//---------------------------------------------------------------------------//
void AdditiveSchwarzWrapper::applyImpl(const MV &x, MV &y) const
{
    // Apply preconditioner
    d_schwarz->apply(x,y,Teuchos::NO_TRANS,1.0,0.0);

    b_num_iters = 1;
}

} // namespace alea

