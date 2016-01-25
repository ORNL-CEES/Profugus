
#ifndef mc_solver_EigenMC_Data_hh
#define mc_solver_EigenMC_Data_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AleaTypedefs.hh"

#include "harness/DBC.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class EigenMC_Data_View
 * \brief Monte Carlo data reorganized into Kokkos Views
 *
 * This class stores matrices needed for MC transport in form of Kokkos::View
 * objects to facilitate device-based parallelism.
 */
//---------------------------------------------------------------------------//
struct EigenMC_Data_View
{
    EigenMC_Data_View(){}

    EigenMC_Data_View( const_scalar_view a,
                  const_scalar_view p,
                  const_scalar_view w,
                  const_ord_view    ind,
                  const_ord_view    off)
        : A(a), P(p), W(w), inds(ind), offsets(off)
    {}

    const_scalar_view A;
    const_scalar_view P;
    const_scalar_view W;
    const_ord_view    inds;
    const_ord_view    offsets;
};

//---------------------------------------------------------------------------//
/*!
 * \class EigenMC_Data_Texture_View
 * \brief Monte Carlo data reorganized into Kokkos texture Views
 *
 * This class stores matrices needed for MC transport in form of Kokkos::View
 * objects to facilitate device-based parallelism.
 */
//---------------------------------------------------------------------------//
struct EigenMC_Data_Texture_View
{
    EigenMC_Data_Texture_View(){}

    EigenMC_Data_Texture_View( EigenMC_Data_View v )
        : A(v.A), P(v.P), W(v.W), inds(v.inds), offsets(v.offsets)
    {}

    random_scalar_view A;
    random_scalar_view P;
    random_scalar_view W;
    random_ord_view    inds;
    random_ord_view    offsets;
};

//---------------------------------------------------------------------------//
/*!
 * \class EigenMC_Data
 * \brief Construct data needed for Monte Carlo linear solver.
 *
 */
//---------------------------------------------------------------------------//

class EigenMC_Data
{
  public:

    // Constructor
    EigenMC_Data(Teuchos::RCP<const MATRIX> A,
            Teuchos::RCP<Teuchos::ParameterList> pl);

    //! Access original matrix, \f$\textbf{A}\f$
    Teuchos::RCP<const MATRIX> getMatrix() const
    {
        REQUIRE( d_Amc != Teuchos::null );
        return d_Amc;
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
    EigenMC_Data_View createKokkosViews();

  private:

    void buildMCMatrix();
    void buildMonteCarloMatrices();
    
    // Original matrix
    Teuchos::RCP<const MATRIX> d_A;

    // matrix manageable for MC
    Teuchos::RCP<CRS_MATRIX> d_Amc;

    // Parameters
    Teuchos::RCP<Teuchos::ParameterList> d_pl;

    // Probability matrix (CDF)
    Teuchos::RCP<CRS_MATRIX> d_P;

    // Weight matrix
    Teuchos::RCP<CRS_MATRIX> d_W;
    
};

}

#endif // mc_solver_EigenMC_Data_hh

