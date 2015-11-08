//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SyntheticAcceleration.hh
 * \author Steven Hamilton
 * \brief  Perform Richardson iteration.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_MultiSplitting_hh
#define mc_solvers_MultiSplitting_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

#include "AleaSolver.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class MultiSplitting
 * \brief Solver linear system using Multi-Splitting techniques
 */
//---------------------------------------------------------------------------//

template <typename solver_type>
struct partition
{
	Teuchos::RCP<solver_type> solver;
	Teuchos::RCP<MV> E;
	Teuchos::ArrayRCP<unsigned int> part;
}


class MultiSplitting : public AleaSolver
{
  public:

    MultiSplitting(Teuchos::RCP<const MATRIX> A,
         Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    // Implementation of apply
    void applyImpl(const MV &x, MV &y) const;

    // Divergence tolerance
    MAGNITUDE d_divergence_tol;
    
    // Parameters for the generations of subdomains
    SCALAR d_num_block;
    SCALAR d_overlap;
    std::string d_inner_solver;
    Teuchos::ArrayRCP< partition<AleaSolver> > subdomains;
};

}

#endif // mc_solvers_SyntheticAcceleration_hh

