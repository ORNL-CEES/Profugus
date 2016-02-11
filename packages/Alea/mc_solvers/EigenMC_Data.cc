//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenMC_Data.cc
 * \author Massimiliano Lupo Pasini
 * \brief  Construct data necessary for Monte Carlo solvers.
 */
//---------------------------------------------------------------------------//

#include <iterator>

#include "EigenMC_Data.hh"
#include "PolynomialBasis.hh"
#include "AleaTypedefs.hh"


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
EigenMC_Data::EigenMC_Data(Teuchos::RCP<const MATRIX> A,
                 Teuchos::RCP<Teuchos::ParameterList> pl)
  : d_A(A)
  , d_pl(pl)
{
    REQUIRE( d_A   != Teuchos::null );
    REQUIRE( d_pl  != Teuchos::null );

    buildMCMatrix();
    buildMonteCarloMatrices();
}

//---------------------------------------------------------------------------//
// Convert matrices to Kokkos View storage
//---------------------------------------------------------------------------//
EigenMC_Data_View EigenMC_Data::createKokkosViews()
{
    REQUIRE( d_Amc->isLocallyIndexed() );
    REQUIRE( d_P->isLocallyIndexed() );
    REQUIRE( d_W->isLocallyIndexed() );
    REQUIRE( d_Amc->supportsRowViews() );
    REQUIRE( d_P->supportsRowViews() );
    REQUIRE( d_W->supportsRowViews() );
    LO numRows     = d_Amc->getNodeNumRows();
    GO numNonZeros = d_Amc->getNodeNumEntries();

    // Allocate views
    scalar_view A("Amc",numNonZeros);
    scalar_view P("Amc",numNonZeros);
    scalar_view W("Amc",numNonZeros);
    ord_view    inds("inds",numNonZeros);
    ord_view    offsets("offsets",numRows+1);

    // Mirror views on host
    scalar_host_mirror A_host       = Kokkos::create_mirror_view(A);
    scalar_host_mirror P_host       = Kokkos::create_mirror_view(P);
    scalar_host_mirror W_host       = Kokkos::create_mirror_view(W);
    ord_host_mirror    inds_host    = Kokkos::create_mirror_view(inds);
    ord_host_mirror    offsets_host = Kokkos::create_mirror_view(offsets);

    Teuchos::ArrayView<const double> A_row, P_row, W_row;
    Teuchos::ArrayView<const int> ind_row;

    // Copy data from CrsMatrix into Kokkos host mirrors
    LO count = 0;
    for( LO irow=0; irow<numRows; ++irow )
    {
        // Get row views
        d_Amc->getLocalRowView(irow,ind_row,A_row);
        d_P->getLocalRowView(irow,ind_row,P_row);
        d_W->getLocalRowView(irow,ind_row,W_row);

        // Copy into host mirror
        std::copy( A_row.begin(), A_row.end(), &A_host(count) );
        std::copy( P_row.begin(), P_row.end(), &P_host(count) );
        std::copy( W_row.begin(), W_row.end(), &W_host(count) );
        std::copy( ind_row.begin(), ind_row.end(), &inds_host(count) );
        count += ind_row.size();
        offsets_host(irow+1) = count;
    }

    // Copy to device
    Kokkos::deep_copy(A,A_host);
    Kokkos::deep_copy(P,P_host);
    Kokkos::deep_copy(W,W_host);
    Kokkos::deep_copy(inds,inds_host);
    Kokkos::deep_copy(offsets,offsets_host);

    // Create data view object
    return EigenMC_Data_View(A,P,W,inds,offsets);
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Build iteration MC matrix Amc = A
//---------------------------------------------------------------------------//
void EigenMC_Data::buildMCMatrix()
{

    // Create Amc
    size_t N = d_A->getNodeNumRows();
    size_t max_nnz = d_A->getNodeMaxNumRowEntries();
    d_Amc = Teuchos::rcp( new CRS_MATRIX(d_A->getDomainMap(),max_nnz) );

    Teuchos::ArrayRCP<SCALAR> Amc_vals(max_nnz), A_vals(max_nnz);
    Teuchos::ArrayRCP<GO>     A_inds(max_nnz);
    size_t num_entries;
    for( size_t irow=0; irow<N; ++irow )
    {
        GO gid = d_A->getDomainMap()->getGlobalElement(irow);

        // Get row from A
        d_A->getGlobalRowCopy(irow,A_inds(),A_vals(),num_entries);

        for( size_t icol=0; icol<num_entries; ++icol )
            Amc_vals[icol] = A_vals[icol];

        /*if(irow < 30) 
	{
	    	std::cout<<irow+1<<"___ ";

		for(int i=0; i<num_entries; ++i)
			std::cout<<A_inds[i]+1<<" ";

		std::cout<<std::endl;
	}*/
        // Add row to Amc
        Teuchos::ArrayView<GO>     ind_view = A_inds(0,num_entries);
        Teuchos::ArrayView<SCALAR> val_view = Amc_vals(0,num_entries);
        d_Amc->insertGlobalValues(irow,ind_view,val_view);
    }

    // Complete construction of Amc
    d_Amc->fillComplete();
    CHECK(d_Amc->isFillComplete());
    CHECK(d_Amc->isStorageOptimized());
}

//---------------------------------------------------------------------------//
// Build probability and weight matrices.
//---------------------------------------------------------------------------//
void EigenMC_Data::buildMonteCarloMatrices()
{
    REQUIRE( d_Amc != Teuchos::null );

    Teuchos::RCP<Teuchos::ParameterList> mc_pl =
        Teuchos::sublist(d_pl,"Monte Carlo");

    // Get absorption probability
    SCALAR abs_prob = mc_pl->get("absorption_probability",0.0);
    SCALAR trans_factor = mc_pl->get("transition_factor",1.0);

    CHECK( d_Amc != Teuchos::null );

    // Now loop over rows in B to build probability/weight matrices
    LO N = d_Amc->getNodeNumRows();
    LO max_nnz = d_Amc->getNodeMaxNumRowEntries();
    d_P = Teuchos::rcp( new CRS_MATRIX(d_Amc->getDomainMap(),max_nnz) );
    d_W = Teuchos::rcp( new CRS_MATRIX(d_Amc->getDomainMap(),max_nnz) );

    Teuchos::ArrayRCP<SCALAR> A_vals(max_nnz), P_vals(max_nnz), W_vals(max_nnz);
    Teuchos::ArrayRCP<GO>     A_inds(max_nnz);
    size_t num_entries;
    for( LO irow=0; irow<N; ++irow )
    {
        d_Amc->getGlobalRowCopy(irow,A_inds(),A_vals(),num_entries);

	//std::cout<<irow<<": ";
        // Compute row sum
        MAGNITUDE row_sum=0.0;
        for( size_t icol=0; icol<num_entries; ++icol )
        {
	    //std::cout<<A_inds[icol]<<" ";
            // Never include zeros
            if( SCALAR_TRAITS::magnitude(A_vals[icol]) > 0.0 )
            {
                row_sum += SCALAR_TRAITS::pow(SCALAR_TRAITS::magnitude(
                            A_vals[icol]),trans_factor);
            }
        }
	//std::cout<<std::endl;
        if( row_sum > 0.0 )
        {
            SCALAR cdf = 0.0;
            for( size_t icol=0; icol<num_entries; ++icol )
            {
                MAGNITUDE this_prob = 0.0;
                if( SCALAR_TRAITS::magnitude(A_vals[icol]) > 0.0 )
                {
                    this_prob = SCALAR_TRAITS::pow(
                        SCALAR_TRAITS::magnitude( A_vals[icol] ), trans_factor )
                        / row_sum;
                }
                this_prob /= (1.0 - abs_prob);
                cdf += this_prob;
                P_vals[icol] = cdf;
                if( this_prob > 0.0 ) 
                    W_vals[icol] = A_vals[icol] / this_prob;
                else
                    W_vals[icol] = 0.0;
            }
        }
        else
        {
            std::fill( P_vals.begin(), P_vals.end(), 0.0 );
            std::fill( W_vals.begin(), W_vals.end(), 0.0 );
        }
        d_P->insertGlobalValues(irow,A_inds(0,num_entries),P_vals(0,num_entries));
        d_W->insertGlobalValues(irow,A_inds(0,num_entries),W_vals(0,num_entries));
    }
    d_P->fillComplete();
    d_W->fillComplete();
}

} // namespace alea

