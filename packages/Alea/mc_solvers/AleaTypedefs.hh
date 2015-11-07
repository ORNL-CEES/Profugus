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
#include "Kokkos_Core.hpp"

namespace alea
{
// These should all be replaced with equivalents in SPn/solvers

//! Type for scalar data
typedef double                                  SCALAR;
//! Type for local indices
typedef int                                     LO;
//! Type for global indices
typedef int                                     GO;
//! Type for Kokkos "Classic" nodes used by Tpetra objects
typedef KokkosClassic::DefaultNode::DefaultNodeType NODE;
//! Device type for Kokkos kernels
typedef Kokkos::DefaultExecutionSpace           DEVICE;
//! Device type for host execution of Kokkos
typedef Kokkos::HostSpace::execution_space      HOST;
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

//type introduced by Max to deal with the partitions of subdomains for MultiSplitting
typedef Teuchos::ArrayRCP< SCALAR >(2, 0) ENDPOINTS;

// Kokkos View types
using Kokkos::View;
typedef Kokkos::MemoryTraits<Kokkos::RandomAccess>        RandomMemory;
typedef DEVICE::scratch_memory_space                      SharedSpace;
typedef Kokkos::MemoryUnmanaged                           UnmanagedMemory;
typedef View<      SCALAR *,DEVICE>                       scalar_view;
typedef View<const SCALAR *,DEVICE>                       const_scalar_view;
typedef View<      LO     *,DEVICE>                       ord_view;
typedef View<const LO     *,DEVICE>                       const_ord_view;
typedef View<const SCALAR *,DEVICE,RandomMemory>          random_scalar_view;
typedef View<const LO     *,DEVICE,RandomMemory>          random_ord_view;
typedef View<      SCALAR *,SharedSpace,UnmanagedMemory>  shared_scalar_view;
typedef View<      LO     *,SharedSpace,UnmanagedMemory>  shared_ord_view;
typedef View<      SCALAR **,DEVICE>                      scalar_view_2d;
typedef View<const SCALAR **,DEVICE>                      const_scalar_view_2d;
typedef View<      SCALAR **,DEVICE,RandomMemory>         random_scalar_view_2d;
typedef View<SCALAR *,DEVICE>::HostMirror                 scalar_host_mirror;
typedef View<LO     *,DEVICE>::HostMirror                 ord_host_mirror;
}

#endif // AleaTypedefs_hh

