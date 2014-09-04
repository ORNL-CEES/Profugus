//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/LinAlgTraits.hh
 * \author Steven Hamilton
 * \brief  Matrices for solver testing.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_test_LinAlgTraits_hh
#define solvers_test_LinAlgTraits_hh

#include <vector>
#include "gtest/utils_gtest.hh"

#include <SPn/config.h>

#include "Teuchos_RCP.hpp"

#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "../TpetraTypedefs.hh"

using profugus::Tpetra_CrsMatrix;
using profugus::Tpetra_MultiVector;

namespace linalg_traits
{

// build_vector
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
    x->PutScalar(0.0);
    return x;
}

template <>
Teuchos::RCP<Tpetra_MultiVector> build_vector<Tpetra_MultiVector>(int N)
{
    return Teuchos::null;
}

// fill_vector
template <class MV>
void fill_vector(Teuchos::RCP<MV> x, std::vector<double> &vals)
{
    NOT_IMPLEMENTED("fill_vector for arbitrary MV type.");
}

template <>
void fill_vector<Epetra_MultiVector>(Teuchos::RCP<Epetra_MultiVector> x,
                                     std::vector<double> &vals)
{
    REQUIRE( vals.size() == x->GlobalLength() );

    for( int i=0; i<x->MyLength(); ++i )
    {
        int global = x->Map().GID(i);
        x->ReplaceGlobalValue(global,0,vals[global]);
    }
}

template <>
void fill_vector<Tpetra_MultiVector>(Teuchos::RCP<Tpetra_MultiVector> x,
                                     std::vector<double> &vals)
{
}

// set_sign
template <class MV>
void set_sign(Teuchos::RCP<MV> x)
{
    NOT_IMPLEMENTED("set_sign for arbitrary MV type.");
}

template <>
void set_sign<Epetra_MultiVector>(Teuchos::RCP<Epetra_MultiVector> x)
{
    double sign = (*x)[0][0] > 0.0 ? 1.0 : -1.0;

    for( int i=0; i<x->MyLength(); ++i )
    {
        (*x)[0][i] = sign * (*x)[0][i];
    }
}

template <>
void set_sign<Tpetra_MultiVector>(Teuchos::RCP<Tpetra_MultiVector> x)
{
}

// test_vector
template <class MV>
void test_vector(Teuchos::RCP<MV> x, std::vector<double> &vals)
{
    NOT_IMPLEMENTED("test_vector for arbitrary MV type.");
}

template <>
void test_vector<Epetra_MultiVector>(Teuchos::RCP<Epetra_MultiVector> x,
                                     std::vector<double> &vals)
{
    REQUIRE( vals.size() == x->GlobalLength() );

    for( int i=0; i<x->MyLength(); ++i )
    {
        int global = x->Map().GID(i);
        EXPECT_SOFTEQ( (*x)[0][i], vals[global], 1e-6 );
    }
}

template <>
void test_vector<Tpetra_MultiVector>(Teuchos::RCP<Tpetra_MultiVector> x,
                                     std::vector<double> &vals)
{
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
    std::vector<std::vector<double> > A_vals =
        {{10.0, 1.1, 2.0, 4.0},
         { 1.1, 9.9, 2.1, 3.2},
         { 0.8, 0.4, 5.3, 1.9},
         { 0.3, 0.1, 0.4, 3.1}};
#ifdef COMM_MPI
    Epetra_MpiComm comm(profugus::communicator);
#else
    Epetra_SerialComm comm;
#endif

    Epetra_Map map( 4, 0, comm );
    Teuchos::RCP<Epetra_CrsMatrix> A( new Epetra_CrsMatrix(Copy,map,4) );
    std::vector<int> inds = {0, 1, 2, 3};

    int local_size = A->NumMyRows();
    int err;
    for( int i=0; i<local_size; ++i )
    {
        int gid = A->GRID(i);
        err = A->InsertGlobalValues(gid,4,&A_vals[gid][0],&inds[0]);
        CHECK( 0 == err );
    }
    A->FillComplete();
    A->OptimizeStorage();
    return A;

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
    std::vector<std::vector<double> > A_vals =
        {{0.56, 0.26, 0.51, 0.26},
         {0.52, 0.13, 0.11, 0.41},
         {0.73, 0.45, 0.40, 0.98},
         {0.30, 0.44, 0.93, 0.35}};
#ifdef COMM_MPI
    Epetra_MpiComm comm(profugus::communicator);
#else
    Epetra_SerialComm comm;
#endif

    Epetra_Map map( 4, 0, comm );
    Teuchos::RCP<Epetra_CrsMatrix> A( new Epetra_CrsMatrix(Copy,map,4) );
    std::vector<int> inds = {0, 1, 2, 3};

    int local_size = A->NumMyRows();
    int err;
    for( int i=0; i<local_size; ++i )
    {
        int gid = A->GRID(i);
        err = A->InsertGlobalValues(gid,4,&A_vals[gid][0],&inds[0]);
        CHECK( 0 == err );
    }
    A->FillComplete();
    A->OptimizeStorage();
    return A;
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
#ifdef COMM_MPI
    Epetra_MpiComm comm(profugus::communicator);
#else
    Epetra_SerialComm comm;
#endif

    Epetra_Map map( N, 0, comm );
    Teuchos::RCP<Epetra_CrsMatrix> A( new Epetra_CrsMatrix(Copy,map,3) );
    std::vector<int> inds(3);
    std::vector<double> vals(3);

    int local_size = A->NumMyRows();
    int err;
    for( int i=0; i<local_size; ++i )
    {
        int gid = A->GRID(i);
        if( gid == 0 )
        {
            inds[0] = 0;
            inds[1] = 1;
            vals[0] = 3.0;
            vals[1] = -1.0;
            err = A->InsertGlobalValues(gid,2,&vals[0],&inds[0]);
            CHECK( 0 == err );
        }
        else if( gid == N-1 )
        {
            inds[0] = N-2;
            inds[1] = N-1;
            vals[0] = -1.0;
            vals[1] = 3.0;
            err = A->InsertGlobalValues(gid,2,&vals[0],&inds[0]);
            CHECK( 0 == err );
        }
        else
        {
            inds[0] = gid-1;
            inds[1] = gid;
            inds[2] = gid+1;
            vals[0] = -1.0;
            vals[1] = 3.0;
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
#ifdef COMM_MPI
    Epetra_MpiComm comm(profugus::communicator);
#else
    Epetra_SerialComm comm;
#endif

    Epetra_Map map( N, 0, comm );
    Teuchos::RCP<Epetra_CrsMatrix> A( new Epetra_CrsMatrix(Copy,map,1) );
    std::vector<int> inds(1);
    std::vector<double> vals(1);

    int local_size = A->NumMyRows();
    int err;
    for( int i=0; i<local_size; ++i )
    {
        int gid = A->GRID(i);
        inds[0] = gid;
        vals[0] = 0.5;
        err = A->InsertGlobalValues(gid,1,&vals[0],&inds[0]);
        CHECK( 0 == err );
    }
    A->FillComplete();
    A->OptimizeStorage();
    return A;
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

} // namespace linalg_traits

#endif // spn_test_LinAlgTraits_hh

//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
