//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/MatrixTraits.hh
 * \author Steven Hamilton
 * \date   Mon Feb 17 13:10:27 2014
 * \brief  MatrixTraits class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_MatrixTraits_hh
#define spn_MatrixTraits_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_OrdinalTraits.hpp"

#include "EpetraExt_RowMatrixOut.h"
#include "MatrixMarket_Tpetra.hpp"

#include "harness/DBC.hh"
#include "solvers/LinAlgTypedefs.hh"


namespace profugus
{

//===========================================================================//
/*!
 * \class MatrixTraits
 * \brief Traits class for Epetra/Tpetra matrices
 */
/*!
 * \example spn/test/tstMatrixTraits.cc
 *
 * Test of MatrixTraits.
 */
//===========================================================================//

template <class T>
class UndefinedMatrixTraits
{
    void NotDefined(){T::this_class_is_missing_a_specialization();}
};

template <class T>
class MatrixTraits
{
  public:
    //@{
    //! Typedefs.
    typedef typename T::MATRIX Matrix_t;
    typedef typename T::MAP    Map_t;
    //@}

    static Teuchos::RCP<Map_t> build_map(int num_local, int num_global,
                                         std::vector<int> indexer )
    {
        UndefinedMatrixTraits<T>::NotDefined();
        return Teuchos::null;
    }

    static Teuchos::RCP<Matrix_t> construct_matrix(
        Teuchos::RCP<const Map_t> map, int num_per_row )
    {
        UndefinedMatrixTraits<T>::NotDefined();
        return Teuchos::null;
    }

    static int local_rows( Teuchos::RCP<const Matrix_t> matrix )
    {
        UndefinedMatrixTraits<T>::NotDefined();
        return Teuchos::null;
    }

    static int global_rows( Teuchos::RCP<const Matrix_t> matrix )
    {
        UndefinedMatrixTraits<T>::NotDefined();
        return Teuchos::null;
    }

    static int global_columns( Teuchos::RCP<const Matrix_t> matrix )
    {
        UndefinedMatrixTraits<T>::NotDefined();
        return Teuchos::null;
    }

    static void add_to_matrix(Teuchos::RCP<Matrix_t> matrix,
                              int row, int count,
                              Teuchos::ArrayRCP<const int> inds,
                              Teuchos::ArrayRCP<const double> vals)
    {
        UndefinedMatrixTraits<T>::NotDefined();
    }

    static void get_local_row_view(Teuchos::RCP<const Matrix_t> matrix, int row,
                                   Teuchos::ArrayView<const int>    &inds,
                                   Teuchos::ArrayView<const double> &vals)
    {
        UndefinedMatrixTraits<T>::NotDefined();
    }

    static int global_col_id(Teuchos::RCP<const Matrix_t> matrix, int local_id)
    {
        UndefinedMatrixTraits<T>::NotDefined();
        return -1;
    }

    static void finalize_matrix(Teuchos::RCP<Matrix_t> matrix)
    {
        UndefinedMatrixTraits<T>::NotDefined();
    }

    static void write_matrix_file(Teuchos::RCP<const Matrix_t> matrix,
                                  std::string filename)
    {
        UndefinedMatrixTraits<T>::NotDefined();
    }

    static UTILS_INT8 global_nonzeros(Teuchos::RCP<const Matrix_t> matrix)
    {
        UndefinedMatrixTraits<T>::NotDefined();
        return -1;
    }
};

// Specialization on EpetraTypes
template <>
class MatrixTraits<EpetraTypes>
{
  public:
    //@{
    //! Typedefs.
    typedef typename EpetraTypes::MATRIX Matrix_t;
    typedef typename EpetraTypes::MAP    Map_t;
    //@}

    static Teuchos::RCP<Map_t> build_map(int num_local, int num_global,
        std::vector<int> indexer=std::vector<int>() )
    {
        // make the communicator
#ifdef COMM_MPI
        Epetra_MpiComm comm(profugus::communicator);
#else
        Epetra_SerialComm comm;
#endif

        // If we have an indexer, use it to construct map, otherwise
        // generate a uniform map.
        Teuchos::RCP<Map_t> map;
        if( !indexer.empty() )
        {
            map = Teuchos::rcp(
                new Map_t(-1, num_local, &indexer[0], 0, comm));
        }
        else
        {
            map = Teuchos::rcp(new Map_t(num_global, 0, comm));
        }

        CHECK(map->NumMyElements() == num_local);
        CHECK(map->NumGlobalElements() == num_global);

        return map;
    }

    static Teuchos::RCP<Matrix_t> construct_matrix(
        Teuchos::RCP<const Map_t> map, int num_per_row )
    {
        Teuchos::RCP<Matrix_t> matrix( new Matrix_t(Copy,*map,num_per_row));
        CHECK(!matrix->StorageOptimized());
        return matrix;
    }

    static int local_rows( Teuchos::RCP<const Matrix_t> matrix )
    {
        return matrix->NumMyRows();
    }

    static int global_rows( Teuchos::RCP<const Matrix_t> matrix )
    {
        return matrix->NumGlobalRows();
    }

    static int global_columns( Teuchos::RCP<const Matrix_t> matrix )
    {
        return matrix->NumGlobalCols();
    }

    static void add_to_matrix(Teuchos::RCP<Matrix_t> matrix,
                              int row, int count,
                              Teuchos::ArrayRCP<const int> inds,
                              Teuchos::ArrayRCP<const double> vals)
    {
        if( count > 0 )
        {
            int err = matrix->InsertGlobalValues(row, count, &vals[0], &inds[0]);
            CHECK( 0 <= err );
        }
    }

    static void get_local_row_view(Teuchos::RCP<const Matrix_t> matrix, int row,
                                   Teuchos::ArrayView<const int>    &inds,
                                   Teuchos::ArrayView<const double> &vals)
    {
        REQUIRE( matrix->IndicesAreLocal() );
        int num_entries;
        int *ind_ptr;
        double *val_ptr;
        int err = matrix->ExtractMyRowView(row,num_entries,val_ptr,ind_ptr);
        CHECK( 0 == err );
        inds = Teuchos::ArrayView<const int>(ind_ptr,num_entries);
        vals = Teuchos::ArrayView<const double>(val_ptr,num_entries);
    }

    static int global_col_id(Teuchos::RCP<const Matrix_t> matrix, int local_id)
    {
        return matrix->GCID(local_id);
    }

    static void finalize_matrix(Teuchos::RCP<Matrix_t> matrix)
    {
        matrix->FillComplete();
        ENSURE(matrix->StorageOptimized());
    }

    static void write_matrix_file(Teuchos::RCP<const Matrix_t> matrix,
                                  std::string filename)
    {
        std::cout << "Writing Epetra matrix to file" << std::endl;
        EpetraExt::RowMatrixToMatrixMarketFile("Epetra.mtx",*matrix,matrix->Label());
    }

    static UTILS_INT8 global_nonzeros(Teuchos::RCP<const Matrix_t> matrix)
    {
        UTILS_INT8 num_nonzeros = matrix->NumMyNonzeros();
        profugus::global_sum( num_nonzeros );
        return num_nonzeros;
    }
};

// Specialization on TpetraTypes
template <>
class MatrixTraits<TpetraTypes>
{
  public:
    //@{
    //! Typedefs.
    typedef typename TpetraTypes::MATRIX Matrix_t;
    typedef typename TpetraTypes::MAP    Map_t;
    //@}

    static Teuchos::RCP<Map_t> build_map(int num_local, int num_global,
        std::vector<int> indexer=std::vector<int>() )
    {
        // make the communicator
        Teuchos::RCP<const Teuchos::Comm<int> > comm =
            Teuchos::DefaultComm<int>::getComm();

        // If we have an indexer, use it to construct map, otherwise
        // generate a uniform map.
        Teuchos::RCP<Map_t> map;
        if( !indexer.empty() )
        {
            map = Teuchos::rcp(
                new Map_t(Teuchos::OrdinalTraits<int>::invalid(),
                          Teuchos::ArrayView<int>(indexer),0,comm) );
        }
        else
        {
            map = Teuchos::rcp(new Map_t(num_global,0,comm) );
        }

        CHECK(map->getNodeNumElements() == num_local);
        CHECK(map->getGlobalNumElements() == num_global);

        return map;
    }

    static Teuchos::RCP<Matrix_t> construct_matrix(
        Teuchos::RCP<const Map_t> map, int num_per_row )
    {
        Teuchos::RCP<Matrix_t> matrix(new Matrix_t(map, num_per_row));
        CHECK(!matrix->isStorageOptimized());
        return matrix;
    }

    static int local_rows( Teuchos::RCP<const Matrix_t> matrix )
    {
        return matrix->getNodeNumRows();
    }

    static int global_rows( Teuchos::RCP<const Matrix_t> matrix )
    {
        return matrix->getGlobalNumRows();
    }

    static int global_columns( Teuchos::RCP<const Matrix_t> matrix )
    {
        return matrix->getGlobalNumCols();
    }

    static void add_to_matrix(Teuchos::RCP<Matrix_t> matrix,
                              int row, int count,
                              Teuchos::ArrayRCP<const int> inds,
                              Teuchos::ArrayRCP<const double> vals)
    {
        matrix->insertGlobalValues(row, inds(0,count), vals(0,count) );
    }

    static void get_local_row_view(Teuchos::RCP<const Matrix_t> matrix, int row,
                                   Teuchos::ArrayView<const int>    &inds,
                                   Teuchos::ArrayView<const double> &vals)
    {
        REQUIRE( matrix->supportsRowViews() );
        REQUIRE( matrix->isLocallyIndexed() );
        matrix->getLocalRowView(row,inds,vals);
    }

    static int global_col_id(Teuchos::RCP<const Matrix_t> matrix, int local_id)
    {
        return matrix->getColMap()->getGlobalElement(local_id);
    }

    static void finalize_matrix(Teuchos::RCP<Matrix_t> matrix)
    {
        matrix->fillComplete();
        ENSURE(matrix->isStorageOptimized());
    }

    static void write_matrix_file(Teuchos::RCP<const Matrix_t> matrix,
                                  std::string filename)
    {
        std::cout << "Writing Tpetra matrix to file" << std::endl;
        Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
            filename,matrix,"tpetra_matrix");
    }

    static UTILS_INT8 global_nonzeros(Teuchos::RCP<const Matrix_t> matrix)
    {
        UTILS_INT8 num_nonzeros = matrix->getNodeNumEntries();
        profugus::global_sum( num_nonzeros );
        return num_nonzeros;
    }
};

} // end namespace profugus

#endif // spn_MatrixTraits_hh

//---------------------------------------------------------------------------//
//                 end of MatrixTraits.hh
//---------------------------------------------------------------------------//
