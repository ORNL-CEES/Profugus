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
#include "Teuchos_DefaultComm.hpp"

#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Operator.h"
#include "../LinAlgTypedefs.hh"

using profugus::EpetraTypes;
using profugus::TpetraTypes;

namespace linalg_traits
{

// build_vector
template <class T>
Teuchos::RCP<typename T::MV> build_vector(int N)
{
    NOT_IMPLEMENTED("build_vector for arbitrary MV type.");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_MultiVector> build_vector<EpetraTypes>(int N)
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
Teuchos::RCP<TpetraTypes::MV> build_vector<TpetraTypes>(int N)
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();

    int index_base = 0;
    Teuchos::RCP<const TpetraTypes::MAP> map(
        new TpetraTypes::MAP(N,index_base,comm) );
    Teuchos::RCP<TpetraTypes::MV> x( new TpetraTypes::MV(map,1) );
    x->putScalar(0.0);
    return x;
}

// fill_vector
template <class T>
void fill_vector(Teuchos::RCP<typename T::MV> x, std::vector<double> &vals)
{
    NOT_IMPLEMENTED("fill_vector for arbitrary MV type.");
}

template <>
void fill_vector<EpetraTypes>(Teuchos::RCP<Epetra_MultiVector> x,
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
void fill_vector<TpetraTypes>(Teuchos::RCP<TpetraTypes::MV> x,
                                     std::vector<double> &vals)
{
    Teuchos::ArrayRCP<double> x_data = x->getDataNonConst(0);
    for( int i=0; i<x->getLocalLength(); ++i )
    {
        int gid = x->getMap()->getGlobalElement(i);
        x_data[i] = vals[gid];
    }
}

// set_sign
template <class T>
void set_sign(Teuchos::RCP<typename T::MV> x)
{
    NOT_IMPLEMENTED("set_sign for arbitrary MV type.");
}

template <>
void set_sign<EpetraTypes>(Teuchos::RCP<Epetra_MultiVector> x)
{
    double sign = (*x)[0][0] > 0.0 ? 1.0 : -1.0;

    for( int i=0; i<x->MyLength(); ++i )
    {
        (*x)[0][i] = sign * (*x)[0][i];
    }
}

template <>
void set_sign<TpetraTypes>(Teuchos::RCP<TpetraTypes::MV> x)
{
    Teuchos::ArrayRCP<double> x_data = x->getDataNonConst(0);
    double sign = x_data[0] > 0.0 ? 1.0 : -1.0;
    for( int i=0; i<x->getLocalLength(); ++i )
    {
        x_data[i] = sign * x_data[i];
    }
}

// test_vector
template <class T>
void test_vector(Teuchos::RCP<typename T::MV> x, std::vector<double> &vals)
{
    NOT_IMPLEMENTED("test_vector for arbitrary MV type.");
}

template <>
void test_vector<EpetraTypes>(Teuchos::RCP<Epetra_MultiVector> x,
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
void test_vector<TpetraTypes>(Teuchos::RCP<TpetraTypes::MV> x,
                                     std::vector<double> &vals)
{
    Teuchos::ArrayRCP<const double> x_data = x->getData(0);
    for( int i=0; i<x->getLocalLength(); ++i )
    {
        int gid = x->getMap()->getGlobalElement(i);
        EXPECT_SOFTEQ( x_data[i], vals[gid], 1e-6 );
    }
}

// build_laplacian
template <class T>
Teuchos::RCP<typename T::MATRIX> build_laplacian(int N)
{
    NOT_IMPLEMENTED("build_laplacian for arbitrary matrix type");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_CrsMatrix> build_laplacian<EpetraTypes>(int N)
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
Teuchos::RCP<TpetraTypes::MATRIX> build_laplacian<TpetraTypes>(int N)
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();

    int index_base = 0;
    Teuchos::RCP<const TpetraTypes::MAP> map(
        new TpetraTypes::MAP(N,index_base,comm) );

    Teuchos::RCP<TpetraTypes::MATRIX> A = Tpetra::createCrsMatrix<double>(map,3);
    Teuchos::ArrayRCP<double> vals(3);
    Teuchos::ArrayRCP<int> inds(3);
    for( int i=0; i<map->getNodeNumElements(); ++i )
    {
        int gid = map->getGlobalElement(i);
        if( 0 == gid )
        {
            inds[0] = 0;
            inds[1] = 1;
            vals[0] = 2.0;
            vals[1] = -1.0;
            A->insertGlobalValues(gid,inds(0,2),vals(0,2));
        }
        else if( N-1 == gid )
        {
            inds[0] = N-2;
            inds[1] = N-1;
            vals[0] = -1.0;
            vals[1] = 2.0;
            A->insertGlobalValues(gid,inds(0,2),vals(0,2));
        }
        else
        {
            inds[0] = gid-1;
            inds[1] = gid;
            inds[2] = gid+1;
            vals[0] = -1.0;
            vals[1] = 2.0;
            vals[2] = -1.0;
            A->insertGlobalValues(gid,inds(),vals());
        }
    }
    A->fillComplete();
    return A;
}

// build_diagonal
template <class T>
Teuchos::RCP<typename T::MATRIX> build_diagonal(int N)
{
    NOT_IMPLEMENTED("build_diagonal for arbitrary matrix type");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_CrsMatrix> build_diagonal<EpetraTypes>(int N)
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
Teuchos::RCP<TpetraTypes::MATRIX> build_diagonal<TpetraTypes>(int N)
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();

    int index_base = 0;
    Teuchos::RCP<const TpetraTypes::MAP> map(
        new TpetraTypes::MAP(N,index_base,comm) );

    Teuchos::RCP<TpetraTypes::MATRIX> A = Tpetra::createCrsMatrix<double>(map,1);
    Teuchos::ArrayRCP<double> vals(1);
    Teuchos::ArrayRCP<int> inds(1);
    for( int i=0; i<map->getNodeNumElements(); ++i )
    {
        int gid = map->getGlobalElement(i);
        inds[0] = gid;
        vals[0] = static_cast<double>(gid+1);
        A->insertGlobalValues(gid,inds(),vals());
    }
    A->fillComplete();
    return A;
}

// build_4x4_lhs
template <class T>
Teuchos::RCP<typename T::MATRIX> build_4x4_lhs()
{
    NOT_IMPLEMENTED("build_4x4_lhs for arbitrary matrix type");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_CrsMatrix> build_4x4_lhs<EpetraTypes>()
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
Teuchos::RCP<TpetraTypes::MATRIX> build_4x4_lhs<TpetraTypes>()
{
    std::vector<std::vector<double> > A_vals =
        {{10.0, 1.1, 2.0, 4.0},
         { 1.1, 9.9, 2.1, 3.2},
         { 0.8, 0.4, 5.3, 1.9},
         { 0.3, 0.1, 0.4, 3.1}};
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();

    int N = 4;
    int index_base = 0;
    Teuchos::RCP<const TpetraTypes::MAP> map(
        new TpetraTypes::MAP(N,index_base,comm) );

    Teuchos::RCP<TpetraTypes::MATRIX> A = Tpetra::createCrsMatrix<double>(map,4);
    Teuchos::ArrayRCP<double> vals;
    std::vector<int> inds_vec = {0, 1, 2, 3};
    Teuchos::ArrayRCP<int> inds = Teuchos::arcp( Teuchos::rcpFromRef(inds_vec) );
    for( int i=0; i<map->getNodeNumElements(); ++i )
    {
        int gid = map->getGlobalElement(i);
        vals = Teuchos::arcp(Teuchos::rcpFromRef(A_vals[gid]));
        A->insertGlobalValues(gid,inds(),vals());
    }
    A->fillComplete();
    return A;
}

// build_4x4_rhs
template <class T>
Teuchos::RCP<typename T::MATRIX> build_4x4_rhs()
{
    NOT_IMPLEMENTED("build_4x4_rhs for arbitrary matrix type");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_CrsMatrix> build_4x4_rhs<EpetraTypes>()
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
Teuchos::RCP<TpetraTypes::MATRIX> build_4x4_rhs<TpetraTypes>()
{
    std::vector<std::vector<double> > A_vals =
        {{0.56, 0.26, 0.51, 0.26},
         {0.52, 0.13, 0.11, 0.41},
         {0.73, 0.45, 0.40, 0.98},
         {0.30, 0.44, 0.93, 0.35}};
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();

    int N = 4;
    int index_base = 0;
    Teuchos::RCP<const TpetraTypes::MAP> map(
        new TpetraTypes::MAP(N,index_base,comm) );

    Teuchos::RCP<TpetraTypes::MATRIX> A = Tpetra::createCrsMatrix<double>(map,4);
    Teuchos::ArrayRCP<double> vals;
    std::vector<int> inds_vec = {0, 1, 2, 3};
    Teuchos::ArrayRCP<int> inds = Teuchos::arcp( Teuchos::rcpFromRef(inds_vec) );
    for( int i=0; i<map->getNodeNumElements(); ++i )
    {
        int gid = map->getGlobalElement(i);
        vals = Teuchos::arcp(Teuchos::rcpFromRef(A_vals[gid]));
        A->insertGlobalValues(gid,inds(),vals());
    }
    A->fillComplete();
    return A;
}

// build_shifted_laplacian
template <class T>
Teuchos::RCP<typename T::MATRIX> build_shifted_laplacian(int N)
{
    NOT_IMPLEMENTED("build_shifted_laplacian for arbitrary matrix type");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_CrsMatrix> build_shifted_laplacian<EpetraTypes>(int N)
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
Teuchos::RCP<TpetraTypes::MATRIX> build_shifted_laplacian<TpetraTypes>(int N)
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();

    int index_base = 0;
    Teuchos::RCP<const TpetraTypes::MAP> map(
        new TpetraTypes::MAP(N,index_base,comm) );

    Teuchos::RCP<TpetraTypes::MATRIX> A = Tpetra::createCrsMatrix<double>(map,3);
    Teuchos::ArrayRCP<double> vals(3);
    Teuchos::ArrayRCP<int> inds(3);
    for( int i=0; i<map->getNodeNumElements(); ++i )
    {
        int gid = map->getGlobalElement(i);
        if( 0 == gid )
        {
            inds[0] = 0;
            inds[1] = 1;
            vals[0] = 3.0;
            vals[1] = -1.0;
            A->insertGlobalValues(gid,inds(0,2),vals(0,2));
        }
        else if( N-1 == gid )
        {
            inds[0] = N-2;
            inds[1] = N-1;
            vals[0] = -1.0;
            vals[1] = 3.0;
            A->insertGlobalValues(gid,inds(0,2),vals(0,2));
        }
        else
        {
            inds[0] = gid-1;
            inds[1] = gid;
            inds[2] = gid+1;
            vals[0] = -1.0;
            vals[1] = 3.0;
            vals[2] = -1.0;
            A->insertGlobalValues(gid,inds(),vals());
        }
    }
    A->fillComplete();
    return A;
}

// build_scaled_identity
template <class T>
Teuchos::RCP<typename T::MATRIX> build_scaled_identity(int N)
{
    NOT_IMPLEMENTED("build_scaled_identity for arbitrary matrix type");
    return Teuchos::null;
}

template <>
Teuchos::RCP<Epetra_CrsMatrix> build_scaled_identity<EpetraTypes>(int N)
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
Teuchos::RCP<TpetraTypes::MATRIX> build_scaled_identity<TpetraTypes>(int N)
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();

    int index_base = 0;
    Teuchos::RCP<const TpetraTypes::MAP> map(
        new TpetraTypes::MAP(N,index_base,comm) );

    Teuchos::RCP<TpetraTypes::MATRIX> A = Tpetra::createCrsMatrix<double>(map,1);
    Teuchos::ArrayRCP<double> vals(1);
    Teuchos::ArrayRCP<int> inds(1);
    for( int i=0; i<map->getNodeNumElements(); ++i )
    {
        int gid = map->getGlobalElement(i);
        inds[0] = gid;
        vals[0] = 0.5;
        A->insertGlobalValues(gid,inds(),vals());
    }
    A->fillComplete();
    return A;
}

template <class T>
Teuchos::RCP<typename T::MATRIX> build_matrix(std::string mat_name, int N)
{
    if( mat_name == "laplacian" )
    {
        return build_laplacian<T>(N);
    }
    else if( mat_name == "diagonal" )
    {
        return build_diagonal<T>(N);
    }
    else if( mat_name == "4x4_lhs" )
    {
        INSIST( 4 == N, "Matrix only defined for N=4" );
        return build_4x4_lhs<T>();
    }
    else if( mat_name == "4x4_rhs" )
    {
        INSIST( 4 == N, "Matrix only defined for N=4" );
        return build_4x4_rhs<T>();
    }
    else if( mat_name == "shifted_laplacian" )
    {
        return build_shifted_laplacian<T>(N);
    }
    else if( mat_name == "scaled_identity" )
    {
        return build_scaled_identity<T>(N);
    }
    return Teuchos::null;
}

} // namespace linalg_traits

#endif // spn_test_LinAlgTraits_hh

//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
