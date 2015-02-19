//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AleaTypedefs.hh
 * \author Steven Hamilton
 * \brief  Common location for frequent typedefs and constants
 */
//---------------------------------------------------------------------------//

#ifndef AleaTypedefs_hh
#define AleaTypedefs_hh

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Teuchos_ScalarTraits.hpp"

// "Classic" Kokkos nodes
#include "Kokkos_DefaultNode.hpp"

// "New" Kokkos nodes
#include "Kokkos_Serial.hpp"
#include "Kokkos_OpenMP.hpp"
#include "Kokkos_Threads.hpp"

namespace alea
{

//! Type for scalar data
typedef double                                  SCALAR;
//! Type for local indices
typedef int                                     LO;
//! Type for global indices
typedef int                                     GO;
//! Type for Kokkos "Classic" nodes used by Tpetra objects
typedef KokkosClassic::DefaultNode::DefaultNodeType NODE;
//! Device type for Kokkos kernels
typedef Kokkos::Threads                         DEVICE;
//! Map type for data distribution
typedef Tpetra::Map<LO,GO,NODE>                 MAP;
//! Multivector type
typedef Tpetra::MultiVector<SCALAR,LO,GO,NODE>  MV;
//! Vector type
typedef Tpetra::Vector<SCALAR,LO,GO,NODE>       VECTOR;
//! Operator type
typedef Tpetra::Operator<SCALAR,LO,GO,NODE>     OP;
//! Matrix type
typedef Tpetra::RowMatrix<SCALAR,LO,GO,NODE>    MATRIX;
//! Concrete matrix type needed for some interfaces
typedef Tpetra::CrsMatrix<SCALAR,LO,GO,NODE>    CRS_MATRIX;
//! Traits class for specified scalar
typedef Teuchos::ScalarTraits<SCALAR>           SCALAR_TRAITS;
//! Traits class for specified local indices
typedef Teuchos::OrdinalTraits<LO>              LO_TRAITS;
//! Traits class for specified global indices
typedef Teuchos::OrdinalTraits<GO>              GO_TRAITS;
//! Type specifying magnitude corresponding to specified scalar
typedef SCALAR_TRAITS::magnitudeType            MAGNITUDE;

}

#endif // AleaTypedefs_hh

