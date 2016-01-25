//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MonteCarloEigenSolver.hh
 * \author Massimiliano Lupo Pasini
 * \brief  Perform adjoint MC histories to solve linear system.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solver_MonteCarloEigenSolver_hh
#define mc_solver_MonteCarloEigenSolver_hh

#include "Kokkos_View.hpp"
#include "Kokkos_Random.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ArrayView.hpp"

#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "AnasaziTypes.hpp"

#include "AleaSolver.hh"
#include "AleaTypedefs.hh"
#include "EigenMC_Data.hh"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class MonteCarloEigenSolver
 * \brief Monte Carlo based eigenvalue solver.
*/
//---------------------------------------------------------------------------//
class MonteCarloEigenSolver : public AleaSolver,
                              public Ifpack2::Preconditioner<SCALAR,LO,GO,NODE>,
                              public Ifpack2::Details::CanChangeMatrix<MATRIX>
{
  public:

    typedef Kokkos::Random_XorShift64_Pool<DEVICE>  generator_pool;

    MonteCarloEigenSolver(Teuchos::RCP<const MATRIX>           A,
                     Teuchos::RCP<Teuchos::ParameterList> pl);

    virtual ~MonteCarloEigenSolver(){}

    // Ifpack2 interface functions

    void setParameters(const Teuchos::ParameterList &list) override
    {
        b_pl->setParameters(list);
    }
    void initialize() override;
    bool isInitialized() const override {return d_initialized;}
    void compute() override { if( !d_initialized ) initialize(); }
    bool isComputed() const override {return d_initialized;}

    Teuchos::RCP<const MATRIX> getMatrix() const override { return b_A; }
    int getNumInitialize() const override { return d_init_count; }
    int getNumCompute() const override { return d_init_count; }
    int getNumApply() const override { return d_apply_count; }
    double getInitializeTime() const override { return 0.0; }
    double getComputeTime() const override { return 0.0; }
    double getApplyTime() const override { return 0.0; }

    // Ifpack2 CanChangeMatrix mixin interface
    void setMatrix(const Teuchos::RCP<const MATRIX> &A) override
    {
        // A can be null, don't check until initialization
        b_A = A;
    }

    // Tpetra::Operator functions that must be redefined due to diamond
    // inheritance

    //! Return domain map of solver (this is range map of operator.)
    virtual Teuchos::RCP<const MAP> getDomainMap() const override
    {
        return b_A->getRangeMap();
    }
    //! Return range map of solver (this is domain map of operator)
    virtual Teuchos::RCP<const MAP> getRangeMap() const override
    {
        return b_A->getDomainMap();
    }

    //! Is transpose apply available.
    virtual bool hasTransposeApply() const override { return false; }

    // Apply
    virtual void apply(const MV &x, MV &y,
                       Teuchos::ETransp mode=Teuchos::NO_TRANS,
                       SCALAR alpha=SCALAR_TRAITS::one(),
                       SCALAR beta=SCALAR_TRAITS::zero()) const override
    {
        // Dispatch to AleaSolver::apply to handle mode,alpha,beta
        AleaSolver::apply(x,y,mode,alpha,beta);
    }

  private:

    // Implementation of apply
    void applyImpl(const MV &x, MV &y) const;

    Teuchos::RCP<Teuchos::ParameterList> d_mc_pl;

    Teuchos::RCP<EigenMC_Data> 	d_mc_data;
    EigenMC_Data_View          	d_mc_data_kokkos;

    enum KERNEL_TYPE { ADAPTIVE, CUDA };

    KERNEL_TYPE d_kernel_type;
    bool    d_use_expected_value;
    GO      d_num_histories;
    SCALAR  d_weight_cutoff;
    SCALAR  d_start_wt_factor;

    // Kokkos random generator pool
    generator_pool d_rand_pool;

    int d_num_threads;

    int  d_init_count;
    bool d_initialized;
    mutable int d_apply_count;
};

} // namespace alea

#endif // mc_solver_MonteCarloEigenSolver_hh

