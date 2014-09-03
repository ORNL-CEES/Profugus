//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/TpetraTypedefs.hh
 * \author Steven Hamilton
 * \date   Mon Jul 01 11:48:09 2013
 * \brief  Collection of typedefs for Tpetra classes.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_TpetraTypedefs_hh
#define solvers_TpetraTypedefs_hh

#include "Teuchos_RCP.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Kokkos_DefaultNode.hpp"

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_Matrix.hpp"

namespace profugus
{

typedef double                                      SCALAR;
typedef int                                         LO;
typedef int                                         GO;
typedef KokkosClassic::DefaultNode::DefaultNodeType NODE;
typedef Tpetra::MultiVector<SCALAR,LO,GO,NODE>      Tpetra_MultiVector;
typedef Tpetra::Vector<SCALAR,LO,GO,NODE>           Tpetra_Vector;
typedef Tpetra::Operator<SCALAR,LO,GO,NODE>         Tpetra_Operator;
typedef Tpetra::RowMatrix<SCALAR,LO,GO,NODE>        Tpetra_RowMatrix;
typedef Tpetra::CrsMatrix<SCALAR,LO,GO,NODE>        Tpetra_CrsMatrix;
typedef Tpetra::Map<LO,GO,NODE>                     Tpetra_Map;
typedef Tpetra::RowGraph<LO,GO,NODE>                Tpetra_Graph;
typedef Xpetra::MultiVector<SCALAR,LO,GO,NODE>      Xpetra_MultiVector;
typedef Xpetra::CrsMatrix<SCALAR,LO,GO,NODE>        Xpetra_CrsMatrix;
typedef Xpetra::Matrix<SCALAR,LO,GO,NODE>           Xpetra_Matrix;

} // end namespace profugus

#endif // solvers_TpetraTypedefs_hh

//---------------------------------------------------------------------------//
//              end of solvers/TpetraTypedefs.hh
//---------------------------------------------------------------------------//
