//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/ChebyshevIteration.hh
 * \author Steven Hamilton
 * \brief  Perform Chebyshev iteration.
 */
//---------------------------------------------------------------------------//

#ifndef Alea_mc_solvers_ChebyshevIteration_hh
#define Alea_mc_solvers_ChebyshevIteration_hh

#include "AleaSolver.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class ChebyshevIteration
 * \brief Perform Chebyshev iteration to solve a linear system
 *
 * The iteration in this class is the implementation given by Algorithm 4
 * of Gutknecht & Roellin, "The Chebyshev iteration revisited."
 * Specifically, updates to the solution vector are made using a "delta"
 * formulation and residuals are computed explicitly.
 */
//---------------------------------------------------------------------------//

class ChebyshevIteration : public AleaSolver
{
  public:

    ChebyshevIteration(Teuchos::RCP<const MATRIX> A,
                        Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    // Implementation of apply
    void applyImpl(const MV &x, MV &y) const;

    SCALAR d_c;
    SCALAR d_d;
};

}

#endif // Alea_mc_solvers_ChebyshevIteration_hh

