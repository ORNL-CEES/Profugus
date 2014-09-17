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

    static void add_to_matrix(Teuchos::RCP<Matrix_t> matrix,
                              int row, int count,
                              Teuchos::ArrayRCP<const int> inds,
                              Teuchos::ArrayRCP<const double> vals)
    {
        UndefinedMatrixTraits<T>::NotDefined();
    }

    static void finalize_matrix(Teuchos::RCP<Matrix_t> matrix)
    {
        UndefinedMatrixTraits<T>::NotDefined();
    }

    static UTILS_INT8 global_nonzeros(Teuchos::RCP<Matrix_t> matrix)
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

    static void add_to_matrix(Teuchos::RCP<Matrix_t> matrix,
                              int row, int count,
                              Teuchos::ArrayRCP<const int> inds,
                              Teuchos::ArrayRCP<const double> vals)
    {
        if( count > 0 )
            matrix->InsertGlobalValues(row, count, &vals[0], &inds[0]);
    }

    static void finalize_matrix(Teuchos::RCP<Matrix_t> matrix)
    {
        matrix->FillComplete();
        ENSURE(matrix->StorageOptimized());
    }

    static UTILS_INT8 global_nonzeros(Teuchos::RCP<Matrix_t> matrix)
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
            map = Teuchos::rcp(
                new Map_t(Teuchos::OrdinalTraits<int>::invalid(),
                          num_global,0,comm) );
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

    static void add_to_matrix(Teuchos::RCP<Matrix_t> matrix,
                              int row, int count,
                              Teuchos::ArrayRCP<const int> inds,
                              Teuchos::ArrayRCP<const double> vals)
    {
        matrix->insertGlobalValues(row, inds(0,count), vals(0,count) );
    }

    static void finalize_matrix(Teuchos::RCP<Matrix_t> matrix)
    {
        matrix->fillComplete();
        ENSURE(matrix->isStorageOptimized());
    }

    static UTILS_INT8 global_nonzeros(Teuchos::RCP<Matrix_t> matrix)
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
