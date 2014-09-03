//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/Test_Matrices.hh
 * \author Steven Hamilton
 * \brief  Matrices for solver testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_test_Test_Matrices_hh
#define solvers_test_Test_Matrices_hh

#include <vector>

#include <SPn/config.h>

#include "Teuchos_RCP.hpp"

#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "../TpetraTypedefs.hh"

using profugus::Tpetra_CrsMatrix;
using profugus::Tpetra_MultiVector;

namespace test_matrix
{

template <class MV>
Teuchos::RCP<MV> build_vector(int N)
{
    NOT_IMPLEMENTED("build_vector for arbitrary MV type.");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_MultiVector> build_vector<Epetra_MultiVector>(int N)
{
#ifdef COMM_MPI
    Epetra_MpiComm comm(profugus::communicator);
#else
    Epetra_SerialComm comm;
#endif

    Epetra_Map map( N, 0, comm );
    Teuchos::RCP<Epetra_MultiVector> x =
        Teuchos::rcp( new Epetra_MultiVector(map,1) );
    x->Random();
    return x;
}

// build_laplacian
template <class Matrix>
Teuchos::RCP<Matrix> build_laplacian(int N)
{
    NOT_IMPLEMENTED("build_laplacian for arbitrary matrix type");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_CrsMatrix> build_laplacian<Epetra_CrsMatrix>(int N)
{
#ifdef COMM_MPI
    Epetra_MpiComm comm(profugus::communicator);
#else
    Epetra_SerialComm comm;
#endif

    Epetra_Map map( N, 0, comm );
    Teuchos::RCP<Epetra_CrsMatrix> A( new Epetra_CrsMatrix(Copy,map,3) );
    std::vector<int> inds(3);
    std::vector<double> vals(3);

    // Boundaries first
    int local_size = A->NumMyRows();
    int err;
    for( int i=0; i<local_size; ++i )
    {
        int gid = A->GRID(i);
        if( gid == 0 )
        {
            inds[0] = 0;
            inds[1] = 1;
            vals[0] = 2.0;
            vals[1] = -1.0;
            err = A->InsertGlobalValues(gid,2,&vals[0],&inds[0]);
            CHECK( 0 == err );
        }
        else if( gid == N-1 )
        {
            inds[0] = N-2;
            inds[1] = N-1;
            vals[0] = -1.0;
            vals[1] = 2.0;
            err = A->InsertGlobalValues(gid,2,&vals[0],&inds[0]);
            CHECK( 0 == err );
        }
        else
        {
            inds[0] = gid-1;
            inds[1] = gid;
            inds[2] = gid+1;
            vals[0] = -1.0;
            vals[1] = 2.0;
            vals[2] = -1.0;
            err = A->InsertGlobalValues(gid,3,&vals[0],&inds[0]);
            CHECK( 0 == err );
        }
    }
    A->FillComplete();
    A->OptimizeStorage();
    return A;
}

template <>
Teuchos::RCP<Tpetra_CrsMatrix> build_laplacian<Tpetra_CrsMatrix>(int N)
{
    return Teuchos::null;
}

// build_diagonal
template <class Matrix>
Teuchos::RCP<Matrix> build_diagonal(int N)
{
    NOT_IMPLEMENTED("build_diagonal for arbitrary matrix type");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_CrsMatrix> build_diagonal<Epetra_CrsMatrix>(int N)
{
#ifdef COMM_MPI
    Epetra_MpiComm comm(profugus::communicator);
#else
    Epetra_SerialComm comm;
#endif

    Epetra_Map map( N, 0, comm );
    Teuchos::RCP<Epetra_CrsMatrix> A( new Epetra_CrsMatrix(Copy,map,1) );
    std::vector<int> inds(1);
    std::vector<double> vals(1);

    // Boundaries first
    int local_size = A->NumMyRows();
    int err;
    for( int i=0; i<local_size; ++i )
    {
        int gid = A->GRID(i);
        inds[0] = gid;
        vals[0] = static_cast<double>(gid+1);
        err = A->InsertGlobalValues(gid,1,&vals[0],&inds[0]);
        CHECK( 0 == err );
    }
    A->FillComplete();
    A->OptimizeStorage();
    return A;
}

template <>
Teuchos::RCP<Tpetra_CrsMatrix> build_diagonal<Tpetra_CrsMatrix>(int N)
{
    return Teuchos::null;
}

// build_4x4_lhs
template <class Matrix>
Teuchos::RCP<Matrix> build_4x4_lhs()
{
    NOT_IMPLEMENTED("build_4x4_lhs for arbitrary matrix type");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_CrsMatrix> build_4x4_lhs<Epetra_CrsMatrix>()
{
    return Teuchos::null;
}

template <>
Teuchos::RCP<Tpetra_CrsMatrix> build_4x4_lhs<Tpetra_CrsMatrix>()
{
    return Teuchos::null;
}

// build_4x4_rhs
template <class Matrix>
Teuchos::RCP<Matrix> build_4x4_rhs()
{
    NOT_IMPLEMENTED("build_4x4_rhs for arbitrary matrix type");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_CrsMatrix> build_4x4_rhs<Epetra_CrsMatrix>()
{
    return Teuchos::null;
}

template <>
Teuchos::RCP<Tpetra_CrsMatrix> build_4x4_rhs<Tpetra_CrsMatrix>()
{
    return Teuchos::null;
}

// build_shifted_laplacian
template <class Matrix>
Teuchos::RCP<Matrix> build_shifted_laplacian(int N)
{
    NOT_IMPLEMENTED("build_shifted_laplacian for arbitrary matrix type");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_CrsMatrix> build_shifted_laplacian<Epetra_CrsMatrix>(int N)
{
    return Teuchos::null;
}

template <>
Teuchos::RCP<Tpetra_CrsMatrix> build_shifted_laplacian<Tpetra_CrsMatrix>(int N)
{
    return Teuchos::null;
}

// build_scaled_identity
template <class Matrix>
Teuchos::RCP<Matrix> build_scaled_identity(int N)
{
    NOT_IMPLEMENTED("build_scaled_identity for arbitrary matrix type");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_CrsMatrix> build_scaled_identity<Epetra_CrsMatrix>(int N)
{
    return Teuchos::null;
}

template <>
Teuchos::RCP<Tpetra_CrsMatrix> build_scaled_identity<Tpetra_CrsMatrix>(int N)
{
    return Teuchos::null;
}

template <class Matrix>
Teuchos::RCP<Matrix> build_matrix(std::string mat_name, int N)
{
    if( mat_name == "laplacian" )
    {
        return build_laplacian<Matrix>(N);
    }
    else if( mat_name == "diagonal" )
    {
        return build_diagonal<Matrix>(N);
    }
    else if( mat_name == "4x4_lhs" )
    {
        INSIST( 4 == N, "Matrix only defined for N=4" );
        return build_4x4_lhs<Matrix>();
    }
    else if( mat_name == "4x4_rhs" )
    {
        INSIST( 4 == N, "Matrix only defined for N=4" );
        return build_4x4_rhs<Matrix>();
    }
    else if( mat_name == "shifted_laplacian" )
    {
        return build_shifted_laplacian<Matrix>(N);
    }
    else if( mat_name == "scaled_identity" )
    {
        return build_scaled_identity<Matrix>(N);
    }
    return Teuchos::null;
}

} // namespace test_matrix

#endif // spn_test_Test_Matrices_hh

//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
