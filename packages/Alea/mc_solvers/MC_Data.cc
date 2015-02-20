//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC_Data.cc
 * \author Steven Hamilton
 * \brief  Construct data necessary for Monte Carlo solvers.
 */
//---------------------------------------------------------------------------//

#include <iterator>

#include "MC_Data.hh"
#include "PolynomialBasis.hh"
#include "AleaTypedefs.hh"

// Trilinos includes
#include "Tpetra_RowMatrixTransposer.hpp"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param A Linear system matrix
 * \param basis PolynomialBasis object defining parameters \f$\alpha,\beta\f$
 * \param pl PL containing relevant parameters
 *
 * The following entries on the "Monte Carlo" sublist of pl are accepted:
 *  - absorption_probability(SCALAR) : scaling parameter applied to each row
 *                                     of probability matrix (1.0)
 */
//---------------------------------------------------------------------------//
MC_Data::MC_Data(Teuchos::RCP<const MATRIX> A,
                 Teuchos::RCP<const PolynomialBasis> basis,
                 Teuchos::RCP<Teuchos::ParameterList> pl)
  : d_A(A)
  , d_pl(pl)
{
    TEUCHOS_ASSERT( d_A   != Teuchos::null );
    TEUCHOS_ASSERT( basis != Teuchos::null );
    TEUCHOS_ASSERT( d_pl  != Teuchos::null );

    buildIterationMatrix(basis);
    buildMonteCarloMatrices();
}


//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Build iteration matrix H=alpha*I + beta*A
//---------------------------------------------------------------------------//
void MC_Data::buildIterationMatrix(Teuchos::RCP<const PolynomialBasis> basis)
{
    TEUCHOS_ASSERT( basis != Teuchos::null );

    // Get parameters for conversion from A to H
    SCALAR alpha, beta;
    basis->getBasisCoefficients(alpha,beta);
    TEUCHOS_ASSERT( alpha != SCALAR_TRAITS::nan() );
    TEUCHOS_ASSERT( beta  != SCALAR_TRAITS::nan() );

    // Create H
    size_t N = d_A->getNodeNumRows();
    size_t max_nnz = d_A->getNodeMaxNumRowEntries();
    d_H = Teuchos::rcp( new CRS_MATRIX(d_A->getDomainMap(),max_nnz) );

    Teuchos::ArrayRCP<SCALAR> H_vals(max_nnz), A_vals(max_nnz);
    Teuchos::ArrayRCP<GO>     A_inds(max_nnz);
    size_t num_entries;
    for( size_t irow=0; irow<N; ++irow )
    {
        GO gid = d_A->getDomainMap()->getGlobalElement(irow);

        // Get row from A
        d_A->getGlobalRowCopy(irow,A_inds(),A_vals(),num_entries);

        for( size_t icol=0; icol<num_entries; ++icol )
        {
            if( A_inds[icol] == gid )
            {
                H_vals[icol] = alpha + beta*A_vals[icol];
            }
            else
            {
                H_vals[icol] = beta*A_vals[icol];
            }
        }

        // Add row to H
        Teuchos::ArrayView<GO>     ind_view = A_inds(0,num_entries);
        Teuchos::ArrayView<SCALAR> val_view = H_vals(0,num_entries);
        d_H->insertGlobalValues(irow,ind_view,val_view);
    }

    // Complete construction of H
    d_H->fillComplete();
    TEUCHOS_ASSERT(d_H->isFillComplete());
    TEUCHOS_ASSERT(d_H->isStorageOptimized());
}

//---------------------------------------------------------------------------//
// Build probability and weight matrices.
//---------------------------------------------------------------------------//
void MC_Data::buildMonteCarloMatrices()
{
    TEUCHOS_ASSERT( d_H != Teuchos::null );

    Teuchos::RCP<Teuchos::ParameterList> mc_pl =
        Teuchos::sublist(d_pl,"Monte Carlo");

    // Get absorption probability
    SCALAR abs_prob = mc_pl->get("absorption_probability",0.0);
    SCALAR trans_factor = mc_pl->get("transition_factor",1.0);

    // Determine if we want forward or adjoint MC
    TEUCHOS_ASSERT(mc_pl->isType<std::string>("mc_type"));
    std::string type = mc_pl->get<std::string>("mc_type");
    TEUCHOS_ASSERT( type == "forward" || type == "adjoint" );

    // Determine "Base" matrix B, such that B = P \circ W
    // Forward MC -> B=H
    // Adjoint MC -> B=H^T
    if( type == "adjoint" )
    {
        // Transpose H
        Tpetra::RowMatrixTransposer<SCALAR,LO,GO,NODE> transposer(d_H);
        d_H = transposer.createTranspose();
    }
    TEUCHOS_ASSERT( d_H != Teuchos::null );

    // Now loop over rows in B to build probability/weight matrices
    LO N = d_H->getNodeNumRows();
    LO max_nnz = d_H->getNodeMaxNumRowEntries();
    d_P = Teuchos::rcp( new CRS_MATRIX(d_H->getDomainMap(),max_nnz) );
    d_W = Teuchos::rcp( new CRS_MATRIX(d_H->getDomainMap(),max_nnz) );

    Teuchos::ArrayRCP<SCALAR> H_vals(max_nnz), P_vals(max_nnz), W_vals(max_nnz);
    Teuchos::ArrayRCP<GO>     H_inds(max_nnz);
    size_t num_entries;
    for( LO irow=0; irow<N; ++irow )
    {
        // Get row from H^T
        d_H->getGlobalRowCopy(irow,H_inds(),H_vals(),num_entries);

        // Compute row sum
        MAGNITUDE row_sum=0.0;
        for( size_t icol=0; icol<num_entries; ++icol )
        {
            // Never include zeros
            if( SCALAR_TRAITS::magnitude(H_vals[icol]) > 0.0 )
            {
                row_sum += SCALAR_TRAITS::pow(SCALAR_TRAITS::magnitude(
                            H_vals[icol]),trans_factor);
            }
        }

        if( row_sum > 0.0 )
        {
            SCALAR cdf = 0.0;
            for( size_t icol=0; icol<num_entries; ++icol )
            {
                MAGNITUDE this_prob = 0.0;
                if( SCALAR_TRAITS::magnitude(H_vals[icol]) > 0.0 )
                {
                    this_prob = SCALAR_TRAITS::pow(
                        SCALAR_TRAITS::magnitude( H_vals[icol] ), trans_factor )
                        / row_sum;
                }
                this_prob /= (1.0 - abs_prob);
                cdf += this_prob;
                P_vals[icol] = cdf;
                if( this_prob > 0.0 )
                {
                    W_vals[icol] = H_vals[icol] / this_prob;
                }
                else
                {
                    W_vals[icol] = 0.0;
                }
            }
        }
        else
        {
            std::fill( P_vals.begin(), P_vals.end(), 0.0 );
            std::fill( W_vals.begin(), W_vals.end(), 0.0 );
        }
        d_P->insertGlobalValues(irow,H_inds(0,num_entries),P_vals(0,num_entries));
        d_W->insertGlobalValues(irow,H_inds(0,num_entries),W_vals(0,num_entries));
    }
    d_P->fillComplete();
    d_W->fillComplete();
}

} // namespace alea

