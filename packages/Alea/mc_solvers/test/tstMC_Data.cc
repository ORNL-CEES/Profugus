//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/test/tstMC_Data.cc
 * \author Steven Hamilton
 * \brief  Test of MC_Data class.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include "gtest/utils_gtest.hh"

#include "../LinearSystem.hh"
#include "../LinearSystemFactory.hh"
#include "../MC_Data.hh"
#include "../PolynomialBasis.hh"
#include "../AleaTypedefs.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace alea;

TEST(MC_Data, Basic)
{
    // Set problem parameters
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(pl,"Problem");
    Teuchos::RCP<Teuchos::ParameterList> poly_pl =
        Teuchos::sublist(pl,"Polynomial");
    Teuchos::RCP<Teuchos::ParameterList> mc_pl =
        Teuchos::sublist(pl,"Monte Carlo");

    // Construct Map
    LO N = 50;
    mat_pl->set("matrix_type","laplacian");
    mat_pl->set("matrix_size",N);
    mat_pl->set("scaling_type","diagonal");

    // Test Adjoint MC
    mc_pl->set("mc_type","adjoint");

    Teuchos::RCP<alea::LinearSystem> system =
        alea::LinearSystemFactory::buildLinearSystem(pl);
    Teuchos::RCP<const MATRIX> A = system->getMatrix();

    Teuchos::RCP<PolynomialBasis> basis( new PolynomialBasis("neumann") );

    Teuchos::RCP<MC_Data> data(new MC_Data(A,basis,pl));
    EXPECT_TRUE( data != Teuchos::null );

    Teuchos::RCP<const MATRIX> H = data->getIterationMatrix();
    EXPECT_TRUE( H != Teuchos::null );
    EXPECT_EQ(H->getGlobalNumRows(), static_cast<size_t>(N));
    EXPECT_EQ(H->getGlobalNumCols(), static_cast<size_t>(N));

    Teuchos::RCP<const MATRIX> P = data->getProbabilityMatrix();
    EXPECT_TRUE( P != Teuchos::null );
    EXPECT_EQ(P->getGlobalNumRows(), static_cast<size_t>(N));
    EXPECT_EQ(P->getGlobalNumCols(), static_cast<size_t>(N));

    Teuchos::RCP<const MATRIX> W = data->getWeightMatrix();
    EXPECT_TRUE( W != Teuchos::null );
    EXPECT_EQ(W->getGlobalNumRows(), static_cast<size_t>(N));
    EXPECT_EQ(W->getGlobalNumCols(), static_cast<size_t>(N));

    // Test values in H, P, W
    // We intentially don't test details of how the matrix is stored, but
    //  rather only the action of the matrix on various vectors
    MV x(A->getDomainMap(),1);
    MV y(A->getDomainMap(),1);
    for( LO icol=0; icol<N; ++icol )
    {
        // Set vector to irow column of identity
        x.putScalar(0.0);
        x.replaceLocalValue(icol,0,1.0);

        // Multiply
        H->apply(x,y);

        // Check result
        for( LO irow=0; irow<N; ++irow )
        {
            // Diagonal entries are 0
            if( irow==icol )
            {
                EXPECT_DOUBLE_EQ( 0.0, y.getData(0)[irow] );
            }
            else if( irow==1 && icol==0 )
            {
                EXPECT_DOUBLE_EQ( 1.0/3.0, y.getData(0)[irow] );
            }
            else if( irow==N-2 && icol==N-1 )
            {
                EXPECT_DOUBLE_EQ( 1.0/3.0, y.getData(0)[irow] );
            }
            // First upper and lower diagonal are 0.5
            else if( std::abs(irow-icol) == 1 )
            {
                EXPECT_DOUBLE_EQ( 0.5, y.getData(0)[irow] );
            }
            // Everything else is zero
            else
            {
                EXPECT_DOUBLE_EQ( 0.0, y.getData(0)[irow] );
            }
        }

        // Multiply
        P->apply(x,y);

        // Check result
        for( LO irow=0; irow<N; ++irow )
        {
            if( irow==0 && icol==0 )
            {
                EXPECT_DOUBLE_EQ( 0.0, y.getData(0)[irow] );
            }
            else if( irow==0 && icol==1 )
            {
                EXPECT_DOUBLE_EQ( 1.0, y.getData(0)[irow] );
            }
            else if( irow==1 && icol==0 )
            {
                EXPECT_DOUBLE_EQ( 2.0/5.0, y.getData(0)[irow] );
            }
            else if( irow==1 && icol==1 )
            {
                EXPECT_DOUBLE_EQ( 2.0/5.0, y.getData(0)[irow] );
            }
            else if( irow==1 && icol==2 )
            {
                EXPECT_DOUBLE_EQ( 1.0, y.getData(0)[irow] );
            }
            else if( irow==N-2 && icol==N-3 )
            {
                EXPECT_DOUBLE_EQ( 3.0/5.0, y.getData(0)[irow] );
            }
            else if( irow==N-2 && icol==N-2 )
            {
                EXPECT_DOUBLE_EQ( 3.0/5.0, y.getData(0)[irow] );
            }
            else if( irow==N-2 && icol==N-1 )
            {
                EXPECT_DOUBLE_EQ( 1.0, y.getData(0)[irow] );
            }
            else if( irow==N-1 && icol==N-2 )
            {
                EXPECT_DOUBLE_EQ( 1.0, y.getData(0)[irow] );
            }
            else if( irow==N-1 && icol==N-1 )
            {
                EXPECT_DOUBLE_EQ( 1.0, y.getData(0)[irow] );
            }
            else if( icol - irow == 1 )
            {
                EXPECT_DOUBLE_EQ( 1.0, y.getData(0)[irow] );
            }
            else if( irow - icol == 1 )
            {
                EXPECT_DOUBLE_EQ( 0.5, y.getData(0)[irow] );
            }
            else if( icol == irow )
            {
                EXPECT_DOUBLE_EQ( 0.5, y.getData(0)[irow] );
            }
            else
            {
                EXPECT_DOUBLE_EQ( 0.0, y.getData(0)[irow] );
            }
        }

        // Multiply
        W->apply(x,y);

        // Check result
        for( LO irow=0; irow<N; ++irow )
        {
            if( irow==0 && icol==1 )
            {
                EXPECT_DOUBLE_EQ( 0.5, y.getData(0)[irow] );
            }
            else if( irow==1 && icol==0 )
            {
                EXPECT_DOUBLE_EQ( 5.0/6.0, y.getData(0)[irow] );
            }
            else if( irow==1 && icol==2 )
            {
                EXPECT_DOUBLE_EQ( 5.0/6.0, y.getData(0)[irow] );
            }
            else if( irow==N-2 && icol==N-3 )
            {
                EXPECT_DOUBLE_EQ( 5.0/6.0, y.getData(0)[irow] );
            }
            else if( irow==N-2 && icol==N-1 )
            {
                EXPECT_DOUBLE_EQ( 5.0/6.0, y.getData(0)[irow] );
            }
            else if( irow==N-1 && icol==N-2 )
            {
                EXPECT_DOUBLE_EQ( 0.5, y.getData(0)[irow] );
            }
            else if( std::abs(icol - irow) == 1 )
            {
                EXPECT_DOUBLE_EQ( 1.0, y.getData(0)[irow] );
            }
            else
            {
                EXPECT_DOUBLE_EQ( 0.0, y.getData(0)[irow] );
            }
        }
    }

}

