//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   testLinearSystemFactory.cc
 * \author Steven Hamilton
 * \brief  Test of LinearSystemFactory class.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "../LinearSystem.hh"
#include "../LinearSystemFactory.hh"
#include "../AleaTypedefs.hh"

using namespace alea;

TEST(LinearSystemFactory, Basic)
{
    // Set problem parameters
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(pl,"Problem");

    // Construct Map
    LO N = 50;
    mat_pl->set("matrix_type","laplacian");
    mat_pl->set("matrix_size",N);

    // Build Laplacian matrix
    Teuchos::RCP<LinearSystem> system =
        LinearSystemFactory::buildLinearSystem(pl);
    Teuchos::RCP<const MATRIX> A = system->getMatrix();

    EXPECT_TRUE( A != Teuchos::null );
    EXPECT_TRUE( A->isFillComplete() );
    EXPECT_EQ( static_cast<size_t>(N), A->getGlobalNumRows() );
    EXPECT_EQ( static_cast<size_t>(N), A->getGlobalNumCols() );

    Teuchos::RCP<const MAP> map = A->getDomainMap();
    LO myNumRows = map->getNodeNumElements();

    // Test values in A
    // We intentially don't test details of how the matrix is stored, but
    //  rather only the action of the matrix on various vectors
    VECTOR x(A->getDomainMap());
    VECTOR y(A->getDomainMap());
    for( LO gcol=0; gcol<N; ++gcol )
    {
        // Set vector to gcol column of identity
        x.putScalar(0.0);
        if( map->isNodeGlobalElement(gcol) )
        {
            x.replaceGlobalValue(gcol,1.0);
        }

        // Multiply
        A->apply(x,y);

        // Check result
        for( LO lrow=0; lrow<myNumRows; ++lrow )
        {
            GO grow = map->getGlobalElement(lrow);
            // Diagonal entries are 1
            if( grow==gcol )
            {
                EXPECT_DOUBLE_EQ( 1.0, y.getData()[lrow] );
            }
            // (0,1) and (N-1,N-2) entries are -1/3
            else if( grow==0 && gcol==1 )
            {
                EXPECT_DOUBLE_EQ( -1.0/3.0, y.getData()[lrow] );
            }
            else if( grow==N-1 && gcol==N-2 )
            {
                EXPECT_DOUBLE_EQ( -1.0/3.0, y.getData()[lrow] );
            }
            // First upper and lower diagonal are -0.5
            else if( std::abs(grow-gcol) == 1 )
            {
                EXPECT_DOUBLE_EQ( -0.5, y.getData()[lrow] );
            }
            // Everything else is zero
            else
            {
                EXPECT_DOUBLE_EQ( 0.0, y.getData()[lrow] );
            }
        }
    }
}

