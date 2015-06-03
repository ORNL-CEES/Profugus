
#ifndef mc_solver_MC_Data_hh
#define mc_solver_MC_Data_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

#include "PolynomialBasis.hh"
#include "harness/DBC.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class MC_Data_View
 * \brief Monte Carlo data reorganized into Kokkos Views
 *
 * This class stores matrices needed for MC transport in form of Kokkos::View
 * objects to facilitate device-based parallelism.
 */
//---------------------------------------------------------------------------//
struct MC_Data_View
{
    MC_Data_View(){}

    MC_Data_View( const_scalar_view h,
                  const_scalar_view p,
                  const_scalar_view w,
                  const_ord_view    ind,
                  const_ord_view    off)
        : H(h), P(p), W(w), inds(ind), offsets(off)
    {}

    const_scalar_view H;
    const_scalar_view P;
    const_scalar_view W;
    const_ord_view    inds;
    const_ord_view    offsets;
};

//---------------------------------------------------------------------------//
/*!
 * \class MC_Data_Texture_View
 * \brief Monte Carlo data reorganized into Kokkos texture Views
 *
 * This class stores matrices needed for MC transport in form of Kokkos::View
 * objects to facilitate device-based parallelism.
 */
//---------------------------------------------------------------------------//
struct MC_Data_Texture_View
{
    MC_Data_Texture_View(){}

    MC_Data_Texture_View( MC_Data_View v )
        : H(v.H), P(v.P), W(v.W), inds(v.inds), offsets(v.offsets)
    {}

    random_scalar_view H;
    random_scalar_view P;
    random_scalar_view W;
    random_ord_view    inds;
    random_ord_view    offsets;
};

//---------------------------------------------------------------------------//
/*!
 * \class MC_Data
 * \brief Construct data needed for Monte Carlo linear solver.
 *
 * This class constructs the matrices \f$\textbf{P}\f$ and \f$\textbf{W}\f$
 * needed to perform the adjoint Monte Carlo algorithm on a lienar system.
 * The original matrix \f$\textbf{A}\f$ is first transformed to an
 * "iteration matrix" \f$\textbf{H}=\alpha \textbf{I} + \beta \textbf{A}\f$
 * (for the standard process based on the Richardson iteration matrix,
 * \f$\alpha=1\f$ and \f$\beta=-1\f$).  The matrix  \f$\textbf{H}\f$ is then
 * decomposed into an element-wise product of a probability transition
 * matrix and a weight transition matrix, i.e.
 * \f$\textbf{H} = \textbf{P} \circ \textbf{W}\f$.  The matrix
 * \f$\textbf{P}\f$ is typically defined to be a row-stochastic matrix
 * (i.e. its rows have a 1-norm of unity), but artificial absorption can be
 * introduced to modify this.
 */
//---------------------------------------------------------------------------//

class MC_Data
{
  public:

    // Constructor
    MC_Data(Teuchos::RCP<const MATRIX> A,
            Teuchos::RCP<const PolynomialBasis> basis,
            Teuchos::RCP<Teuchos::ParameterList> pl);

    //! Access original matrix, \f$\textbf{A}\f$
    Teuchos::RCP<const MATRIX> getMatrix() const
    {
        REQUIRE( d_A != Teuchos::null );
        return d_A;
    }
    //! Access iteration matrix, \f$\textbf{H}\f$
    Teuchos::RCP<const MATRIX> getIterationMatrix() const
    {
        REQUIRE( d_H != Teuchos::null );
        return d_H;
    }
    //! Access probability transition matrix, \f$\textbf{P}\f$
    Teuchos::RCP<const MATRIX> getProbabilityMatrix() const
    {
        REQUIRE( d_P != Teuchos::null );
        return d_P;
    }
    //! Access weight transition matrix, \f$\textbf{W}\f$
    Teuchos::RCP<const MATRIX> getWeightMatrix() const
    {
        REQUIRE( d_W != Teuchos::null );
        return d_W;
    }

    //! Convert matrices to Kokkos::Views
    MC_Data_View createKokkosViews();

  private:

    void buildIterationMatrix(Teuchos::RCP<const PolynomialBasis> basis);
    void buildMonteCarloMatrices();

    // Original matrix
    Teuchos::RCP<const MATRIX> d_A;

    // Parameters
    Teuchos::RCP<Teuchos::ParameterList> d_pl;

    // Richardson iteration matrix
    Teuchos::RCP<CRS_MATRIX> d_H;

    // Probability matrix (CDF)
    Teuchos::RCP<CRS_MATRIX> d_P;

    // Weight matrix
    Teuchos::RCP<CRS_MATRIX> d_W;
};

}

#endif // mc_solver_MC_Data_hh

